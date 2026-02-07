classdef MagBox < MagSource
    
    % MagBox: class representing a magnetized rectangular prism.
    %
    %   MagBox( vx, vy, vz, M ) is the standard constructor where:
    %
    %   The dimensions of the box are specified by:
    %       vx (1x2 array): range of x positions spanned by the box.
    %       vy (1x2 array): range of y positions spanned by the box.
    %       vz (1x2 array): range of z positions spanned by the box.
    %
    %   The box is assumed to have uniform magnetization:
    %       M (3x1 vector): remanent magnetization (A/m)
    %
    %   Optional additional properties:
    %       chi (scalar): magnetic susceptibility (n/d); if included, the
    %       total magnetization is computed internally as M = Mr + Mi where
    %       Mr is the remanent magnetization specified above and Mi = Ba *
    %       chi/mu0 is the induced magnetization (additionally requires an
    %       ambient field, Ba, to be specified)
    %       Ba (3x1 vector): ambient magnetic field (Teslas)
    %
    %   Derived properties include:
    %       V (scalar): box volume (m^3)
    %       mtot (scalar): total magnetic moment (Am^2)
    %
    %
    %   Disclaimer: This code is provided as-is, has been tested only very
    %   informally, and may not always behave as intended. I find it useful for
    %   my own work, and I hope you will too, but I make no guarantees as to
    %   the accuracy or robustness of this code. This code is also actively
    %   under development and future versions may not be backward compatible.
    %
    %   Doug Hemingway (djheming@berkeley.edu)
    %   University of California, Berkeley
    %   2018-06-07
    %
    
    properties (Access=private)

        int_vx % (1x2 array) range of x positions spanned by the box in its own frame.
        int_vy % (1x2 array) range of y positions spanned by the box in its own frame.
        int_vz % (1x2 array) range of z positions spanned by the box in its own frame.

    end

    properties (Access=public)
        
        Mframe = 'A' % coordinate frame in which M is expressed (analysis coordinate system by default).
        Mr % (3x1 vector) remanent magnetization (A/m), assumed uniform throughout the box.
        chi % (n/d scalar) magnetic susceptibility for induced magnetization.

        % The vx, vy, vz ranges are specified in a coordinate system
        % that may or may not be unique to this source body. If the source
        % coordinate system (S) is distinct from the analysis coordinate
        % system (A) chosen for the problem at hand, the relationship
        % between the two coordinate systems is captured with the following
        % vector and rotation matrix. 
        v_AS = zeros(3,1) % (3x1 vector) position of source coordinate system with respect to the analysis coordinate system
        R_AS = eye(3) % (3x3 matrix) rotation matrix that transforms from the source coordinate system to the analysis coordinate system
        % By default, the two coordinate systems are identical so the
        % vector is zero and the rotation matrix is the identiy matrix.
        % But, in general, the two may be distinct.  
        % As an example, if the two coordinate systems are distinct, then a
        % position, p, can be expressed in either coordinate system and the
        % two vectors would be related by pA = v_AS + R_AS * pS. Or,
        % equivalently, pS = R_SA * ( pA - v_AS ). Recall that the inverse
        % of a rotation matrix is also its transpose, so R_SA = R_AS'.
        
    end
    
    properties (Dependent)

        M % (3x1 vector) magnetization (A/m), assumed uniform throughout the box.
        Mi % (3x1 vector) induced magnetization (A/m), assumed uniform throughout the box.
        vx % (1x2 array) range of x positions spanned by the box in its own frame.
        vy % (1x2 array) range of y positions spanned by the box in its own frame.
        vz % (1x2 array) range of z positions spanned by the box in its own frame.
        vS % (3x8 array) vector positions of each of the 8 vertices in the prism's own frame (source frame).
        wx % prism width in the x-direction in its own frame.
        wy % prism width in the y-direction in its own frame.
        wz % prism width in the z-direction in its own frame.
        w % (1x3 array) summarizing widths in each dimension, in the source frame.
        x0 % prism's center in x-direction in its own frame.
        y0 % prism's center in y-direction in its own frame.
        z0 % prism's center in z-direction in its own frame.
        cenS % (3x1 vector) summarizing prism's centroid in its own frame (x0,y0,z0).
        finitedims % (logical) identifies which dimensions are finite
        finitefaces % (3x2 logical) identifying which faces are finite

        R_SA % the inverse of the R_AS rotation matrix
        MA % (3x1 vector) magnetization (A/m) in the analysis frame.
        vA % (3x8 array) vector positions of each of the 8 vertices in the analysis frame.
        cenA % (3x1 vector) vector position of prism centroid in the analysis frame.
        xrng % (1x2 array) range of x positions spanned by the box in the analysis frame.
        yrng % (1x2 array) range of y positions spanned by the box in the analysis frame.
        zrng % (1x2 array) range of z positions spanned by the box in the analysis frame.
        rng % (3x2 array) summarizing the above ranges.
        
        mtot % total magnetic moment (Am^2)
        Vtot % total volume of the box (m^3)
        
    end
    
    methods
        
        % Constructor.
        function newMagBox = MagBox( varargin )
            if nargin >= 4
                newMagBox.int_vx = sort(varargin{1});
                newMagBox.int_vy = sort(varargin{2});
                newMagBox.int_vz = sort(varargin{3});
                newMagBox.Mr = varargin{4};
                if nargin >= 6
                    newMagBox.v_AS = varargin{5};
                    newMagBox.R_AS = varargin{6};
                    if nargin >= 7
                        newMagBox.Mframe = varargin{7};
                        if nargin >= 8
                            args = BaseTools.argarray2struct( varargin(8:end) );
                            if isfield(args,'chi')
                                newMagBox.chi = args.chi;
                            end
                        end
                    end
                end
            end
        end

        % Setter methods to ensure faces are sorted in correct order.
        function thisMagBox = set.vx( thisMagBox, new_vx )
            thisMagBox.int_vx = sort(new_vx);
        end
        function thisMagBox = set.vy( thisMagBox, new_vy )
            thisMagBox.int_vy = sort(new_vy);
        end
        function thisMagBox = set.vz( thisMagBox, new_vz )
            thisMagBox.int_vz = sort(new_vz);
        end

        % Allow the user to set M directly.
        function thisMagBox = set.M( thisMagBox, new_M )
            % If the user is trying to directly set the total M, we have to
            % assume they want a remanent magnetization only (i.e., does
            % not depend on the ambient field).
            thisMagBox.Mr = new_M;
            thisMagBox.chi = 0; % This will make Mi = 0.
        end

        % Simple getter methods.
        function Mr = get.Mr( thisMagBox )
            % Return remanent magnetization vector. Zero of not set.
            if isempty( thisMagBox.Mr )
                Mr = zeros(3,1);
            else
                Mr = thisMagBox.Mr;
            end
        end
        function Mi = get.Mi( thisMagBox )
            % Return induced magnetization given ambient field and magnetic
            % susceptibility. Return zero if either of those is missing.
            if ~isempty( thisMagBox.chi ) && ~isempty( thisMagBox.Ba )
                Mi = thisMagBox.Ba * thisMagBox.chi / MagSource.mu_naught;
            else
                Mi = zeros(3,1);
            end
        end
        function M = get.M( thisMagBox )
            % Total magnetization vector is the remanent plus induced M.
            M = thisMagBox.Mr + thisMagBox.Mi;
        end
        function vx = get.vx( thisMagBox )
            vx = thisMagBox.int_vx;
        end
        function vy = get.vy( thisMagBox )
            vy = thisMagBox.int_vy;
        end
        function vz = get.vz( thisMagBox )
            vz = thisMagBox.int_vz;
        end
        function vS = get.vS( thisMagBox )
            vS = zeros(3,2,2,2);
            for a = 1 : 2
                for b = 1 : 2
                    for c = 1 : 2
                        vS(1,a,b,c) = thisMagBox.vx(a);
                        vS(2,a,b,c) = thisMagBox.vy(b);
                        vS(3,a,b,c) = thisMagBox.vz(c);
                    end
                end
            end
            vS = reshape(vS,3,8);
        end
        function vA = get.vA( thisMagBox )
            if isequal( thisMagBox.R_AS, eye(3) )
                vA = thisMagBox.v_AS + thisMagBox.vS;
            elseif any(isinf(thisMagBox.vS(:)))
                error( 'Cannot rotate infinitely long prism.' );
            else
                vA = thisMagBox.v_AS + thisMagBox.R_AS * thisMagBox.vS;
            end
        end
        function MA = get.MA( thisMagBox )
            MA = thisMagBox.R_AS * thisMagBox.M;
        end
        function cenA = get.cenA( thisMagBox )
            cenA = thisMagBox.v_AS + thisMagBox.R_AS * thisMagBox.cenS;
        end
        function xrng = get.xrng( thisMagBox )
            xrng = [ min(thisMagBox.vA(1,:)) max(thisMagBox.vA(1,:)) ];
        end
        function yrng = get.yrng( thisMagBox )
            yrng = [ min(thisMagBox.vA(2,:)) max(thisMagBox.vA(2,:)) ];
        end
        function zrng = get.zrng( thisMagBox )
            zrng = [ min(thisMagBox.vA(3,:)) max(thisMagBox.vA(3,:)) ];
        end
        function rng = get.rng( thisMagBox )
            rng = [ thisMagBox.xrng; thisMagBox.yrng; thisMagBox.zrng ];
        end
        function wx = get.wx( thisMagBox )
            wx = thisMagBox.vx(2)-thisMagBox.vx(1);
        end
        function wy = get.wy( thisMagBox )
            wy = thisMagBox.vy(2)-thisMagBox.vy(1);
        end
        function wz = get.wz( thisMagBox )
            wz = thisMagBox.vz(2)-thisMagBox.vz(1);
        end
        function w = get.w( thisMagBox )
            w = [ thisMagBox.wx thisMagBox.wy thisMagBox.wz ];
        end
        function finitedims = get.finitedims( thisMagBox )
            finitedims = [ ~isinf(thisMagBox.wx) ~isinf(thisMagBox.wy) ~isinf(thisMagBox.wz) ];
        end
        function finitefaces = get.finitefaces( thisMagBox )
            finitefaces = [
                ~isinf(thisMagBox.vx);
                ~isinf(thisMagBox.vy);
                ~isinf(thisMagBox.vz)
                ];
        end
        function x0 = get.x0( thisMagBox )
            if isinf(thisMagBox.wx)
                x0 = 0;
            else
                x0 = ( thisMagBox.vx(1)+thisMagBox.vx(2) )/2;
            end
        end
        function y0 = get.y0( thisMagBox )
            if isinf(thisMagBox.wy)
                y0 = 0;
            else
                y0 = ( thisMagBox.vy(1)+thisMagBox.vy(2) )/2;
            end
        end
        function z0 = get.z0( thisMagBox )
            if isinf(thisMagBox.wz)
                z0 = 0;
            else
                z0 = ( thisMagBox.vz(1)+thisMagBox.vz(2) )/2;
            end
        end
        function cenS = get.cenS( thisMagBox )
            cenS = [ thisMagBox.x0; thisMagBox.y0; thisMagBox.z0 ];
        end
        function R_SA = get.R_SA( thisMagBox )
            R_SA = thisMagBox.R_AS';
        end
        function Vtot = get.Vtot( thisMagBox )
            Vtot = (max(thisMagBox.vx)-min(thisMagBox.vx)) * ...
                (max(thisMagBox.vy)-min(thisMagBox.vy)) * ...
                (max(thisMagBox.vz)-min(thisMagBox.vz));
        end
        function mtot = get.mtot( thisMagBox )
            mtot = thisMagBox.Vtot * norm( thisMagBox.M );
        end
        
        % Compute fields.
        B = computeBfield( thisMagBox, p, varargin )
        Q = computeQfield( thisMagBox, p, varargin )       
        
        % For 2D prisms (infinitely long in one dimension), we have closed
        % form analytical expressions for the zero crossings.
        function p = computeQiiZeros( thisMagBox, u )
            if sum(thisMagBox.finitedims) == 2
                % Sort and index the dimensions from shortest to longest.
                [ ~, ind ] = sort( thisMagBox.w ); % ind=1 for the shorter dimension, ind=2 for the longer dimension, ind=3 for the infinite.
                v = sqrt( u.^2 + ( thisMagBox.w(ind(2))^2 - thisMagBox.w(ind(1))^2 )/4 );
                % Now we have a vector (u) corresponding to the shorter
                % dimension and another (v) corresponding to the longer
                % dimension. We want the other (infinite) dimension to be
                % all zeros). We also want to NaN-out any points inside the
                % prism (i.e., where the long coordinates are less than the
                % longer width of the prism).
                v_inside = v > thisMagBox.rng(ind(2),1) & v < thisMagBox.rng(ind(2),2);
                v(v_inside) = NaN;
                other = zeros(size(u));
                p_unsorted = [ u NaN u; v NaN -v; other NaN other ]; % NaN separators ensure that these are not treated as one unbroken line.
                p = p_unsorted(ind,:);
            else
                warning( 'Can only compute Qii zeros if the prism has exactly two finite dimensions.' );
            end
        end
        function p = computeQijZeros( thisMagBox, u )
            if sum(thisMagBox.finitedims) == 2
                % Sort and index the dimensions from shortest to longest.
                [ ~, ind ] = sort( thisMagBox.w ); % ind=1 for the shorter dimension, ind=2 for the longer dimension, ind=3 for the infinite.
                v = u;
                other = zeros(size(u));
                p_unsorted = [ u NaN other; other NaN v; other NaN other ];
                % Now find all the points that are inside the prism and set
                % them to NaN.
                inside = p_unsorted(1,:) > thisMagBox.rng(ind(1),1) & p_unsorted(1,:) < thisMagBox.rng(ind(1),2) & ...
                    p_unsorted(2,:) > thisMagBox.rng(ind(2),1) & p_unsorted(2,:) < thisMagBox.rng(ind(2),2);
                p_unsorted(:,inside) = NaN;
                p = p_unsorted(ind,:);
            else
                warning( 'Can only compute Qij zeros if the prism has exactly two finite dimensions.' );
            end
        end
        function [ p_short, p_long ] = computeQijPeaks( thisMagBox, u )
            if sum(thisMagBox.finitedims) == 2
                % Sort and index the dimensions from shortest to longest.
                [ ~, ind ] = sort( thisMagBox.w ); % ind=1 for the shorter dimension, ind=2 for the longer dimension, ind=3 for the infinite.
                wu = thisMagBox.w(ind(1));
                wv = thisMagBox.w(ind(2));
                % Find the peaks in d/dv.
                vInnerRoot = sqrt( wv^4 + wv^2 * wu^2 + 4*wv^2 * u.^2 + wu^4 - 4*wu^2 * u.^2 + 16*u.^4 );
                v = sqrt( (1/12)*(wv^2 - wu^2) - (1/3)*u.^2 + (1/6)*vInnerRoot );
                other = zeros(size(u));
                p_unsorted = [ u NaN u; v NaN -v; other NaN other ]; % NaN separators ensure that these are not treated as one unbroken line.
                p_short = p_unsorted(ind,:);
                % Find the peaks in d/du.
                v = u;
                uInnerRoot = sqrt( wu^4 + wu^2 * wv^2 + 4*wu^2 * v.^2 + wv^4 - 4*wv^2 * v.^2 + 16*v.^4 );
                u = sqrt( (1/12)*(wu^2 - wv^2) - (1/3)*v.^2 + (1/6)*uInnerRoot );
                p_unsorted = [ u NaN -u; v NaN v; other NaN other ]; % NaN separators ensure that these are not treated as one unbroken line.
                p_long = p_unsorted(ind,:);
            else
                warning( 'Can only compute Qij peaks if the prism has exactly two finite dimensions.' );
            end
        end

        % Display functions.
        ah = drawBox( thisMagBox, varargin )
        function fh = showM( thisMagBox, varargin )
            % TO DO: make this nicer!
            [ args, optargs ] = BaseTools.argarray2struct(varargin, { 'Color', 'r', 'Color2', [.6 .6 .6], 'show_comps', false, 'comps', 'xyz', 'view', [] });
            lenM = norm(thisMagBox.M);
            if isfield( args, 'M_length' ) && ~isempty( args.M_length )
                lenM = lenM/args.M_length;
            end
            Mdisp_x = thisMagBox.x0+[ 0 thisMagBox.M(1)/lenM ];
            Mdisp_y = thisMagBox.y0+[ 0 thisMagBox.M(2)/lenM ];
            Mdisp_z = thisMagBox.z0+[ 0 thisMagBox.M(3)/lenM ];
            if isfield( args, 'show_comps' ) && args.show_comps
                headlength = 1/10;
                fi = gobjects(size(args.comps));
                for i = 1 : numel(args.comps)
                    [ ai, fi(i) ] = BaseTools.verify_axes_handle;
                    BaseTools.drawArrow( ai, Mdisp_x*(strcmp(args.comps(i),'x')), Mdisp_y*(strcmp(args.comps(i),'y')), Mdisp_z*(strcmp(args.comps(i),'z')), 'LineWidth', 2.0, 'Color', args.Color, 'headlength', headlength );
                    BaseTools.drawArrow( ai, Mdisp_x, Mdisp_y, Mdisp_z, 'LineWidth', 2.0, 'Color', args.Color2, 'headlength', headlength, 'view', args.view );
                end
                fh = BaseTools.tileFigures( fi' );
            else
                [ ah, fh ] = BaseTools.extractAxesHandle( args );
                BaseTools.drawArrow( ah, Mdisp_x, Mdisp_y, Mdisp_z, 'LineWidth', 2.0, optargs{:} );
                if isfield( args, 'view' ) && ~isempty( args.view )
                    view( ah, args.view );
                    camlight; % Give it some 3D lighting.
                end
            end
        end

        % Make a movie from a set of key frames.
        function previewBmovie( keyBoxes, varargin )

            % Parse input arguments.
            args = BaseTools.argarray2struct( varargin, { 'duration', 3, 'fixed_moment', [], 'b', 12 } );
            num_keys = length(keyBoxes);

            % Do we need to normalize the magnetic moments?
            if isfield( args, 'fixed_moment' ) && ~isempty( args.fixed_moment )
                for t = 1 : num_keys
                    keyBoxes(t).M = keyBoxes(t).M * args.fixed_moment/keyBoxes(t).mtot;
                end
            end
            survey_volume = SurveyField( linspace(-args.b,args.b), linspace(-args.b,args.b), linspace(-args.b,args.b) );

            % If the prism geometry is not changing, we can pre-compute Q
            % and avoid having to re-compute it every time we compute B.
            if consistent_geometry( keyBoxes )
                Q = keyBoxes(1).computeQfield( survey_volume.p );
                fprintf( 'Prism geometry is fixed so we can pre-compute Q.\n' );
            else
                Q = [];
            end

            % Display all of the key frames.
            for t = 1 : num_keys
                myBox = keyBoxes(t);
                BaseTools.tileFigures( myBox.showBfieldContours( 'xyz', survey_volume, 'Q', Q, 'view', [ 25 15 ], 'title', true ), 1, 3 );
            end
            
        end        
        function makeBmovie( keyBoxes, output_filename, varargin )

            % Parse input arguments.
            args = BaseTools.argarray2struct( varargin, { 'duration', 3, 'fixed_moment', [], 'b', 12 } );
            fprintf( 'Generating %.1f-second movie called: ''%s''...\n', args.duration, output_filename );

            % Figure out how many frames we will have and where the key frames should
            % appear in that sequence.
            num_keys = length(keyBoxes);
            if num_keys == 1
                error( 'Need more than one key frame to make a movie.' );
            end
            num_frames = args.duration * 30;
            [ frame_inds, key_frame_inds ] = BaseTools.spaced_sequence( num_keys, num_frames );
            boxList = MagBox.interpBox( key_frame_inds, keyBoxes, frame_inds );

            % Do we need to normalize the magnetic moments?
            if isfield( args, 'fixed_moment' ) && ~isempty( args.fixed_moment )
                for t = 1 : num_frames
                    boxList(t).M = boxList(t).M * args.fixed_moment/boxList(t).mtot;
                end
            end
            survey_volume = SurveyField( linspace(-args.b,args.b), linspace(-args.b,args.b), linspace(-args.b,args.b) );

            % If the prism geometry is not changing, we can pre-compute Q
            % and avoid having to re-compute it every time we compute B.
            if consistent_geometry( keyBoxes )
                Q = keyBoxes(1).computeQfield( survey_volume.p );
                fprintf( 'Prism geometry is fixed so we can pre-compute Q.\n' );
            else
                Q = [];
            end

            % Set up movie file and writer object.
            writerObj = VideoWriter( output_filename, 'MPEG-4' );
            writerObj.Quality = 95;
            open( writerObj );

            % Make a series of frames from the series of box objects.
            movietimer = tic;
            for t = 1 : num_frames
                myBox = boxList(t);
                comp_fig = BaseTools.tileFigures( myBox.showBfieldContours( 'xyz', survey_volume, 'Q', Q, 'view', [ 25 15 ], 'title', true ), 1, 3 );
                writeVideo( writerObj, getframe(comp_fig) );
                close(comp_fig);
                BaseTools.progress_report(movietimer,t,num_frames);
            end

            % Save movie file.
            close( writerObj );
            fprintf( 'Completed movie in %s. Saved ''%s''.\n', BaseTools.timetostr(toc(movietimer)), output_filename );

        end

        % Check whether prisms have consistent geometry.
        function result = consistent_geometry( boxList )
            result = true; % Start by assuming they're all the same and then return false if we find any exceptions.
            N = length(boxList);
            if N > 1
                for k = 2 : N
                    if any( boxList(k).xrng ~= boxList(1).xrng ) || ...
                            any( boxList(k).yrng ~= boxList(1).yrng ) || ...
                            any( boxList(k).zrng ~= boxList(1).zrng )
                        result = false;
                        return;
                    end
                end
            end
        end
        
    end

    methods (Static)

        % Core functions for computing elements of Q.
        [ Qii, alphas, thetas ] = computeQii( ri, rj, rk )
        [ Qij, f_ab, f_c, ts, Rsq ] = computeQij( ri, rj, rk )
        
        % Basic calculations.
        function y_zero = compute_Qii_zeros( x, wx, wy )
            y_zero = sqrt( x.^2 + ( wy^2 - wx^2 )/4 );
        end
        function y_peak = compute_Qij_peaks( x, wx, wy )
            cornersqdiff = wy^2 - wx^2;
            yInnerRoot = sqrt( wy^4 + wy^2 * wx^2 + wx^4 + 4*cornersqdiff*x.^2 + 16*x.^4 );
            y_peak = sqrt( cornersqdiff/12 - (1/3)*x.^2 + (1/6)*yInnerRoot );
        end

        % Interpolate from a list of MagBox objects.
        function outputBoxList = interpBox( x, inputBoxList, xq )
            % The x and xq vectors are arbitrary indeces but xq must lie
            % within the range of x. That is, xq can be any set of real
            % numbers >=x(1) and <=x(end). 
            num_input_boxes = length(inputBoxList);
            num_output_boxes = length(xq);
            if num_input_boxes ~= length(x)
                error( 'x and the input box list must have the same lengths.' );
            end
            if max(xq) > max(x) || min(xq) < min(x)
                error( 'xq must lie within the range of x.' );
            end
            outputBoxList = repmat( MagBox, num_output_boxes, 1 );
            % Construct an array to hold all the data. 
            boxdims = nan(num_input_boxes,6);
            boxMs = nan(num_input_boxes,3);
            for i = 1 : num_input_boxes
                boxdims(i,:) = [ inputBoxList(i).vx inputBoxList(i).vy inputBoxList(i).vz ];
                boxMs(i,:) = inputBoxList(i).M';
            end
            % For the box dimensions, make a piece-wise linear interpolation.
            output_boxdims = interp1( x, boxdims, xq );
            for k = 1 : num_output_boxes
                outputBoxList(k).vx = output_boxdims(k,1:2);
                outputBoxList(k).vy = output_boxdims(k,3:4);
                outputBoxList(k).vz = output_boxdims(k,5:6);
            end
            % For the magnetization vector, do spherical vector interpolation.
            outputM = BaseTools.interp_vectors( x, boxMs', xq );
            for k = 1 : num_output_boxes
                outputBoxList(k).M = outputM(:,k);
            end
        end

        % Testing.
        function varargout = unit_test( type )

            close all;
            if ~exist( 'type', 'var' )
                type = 'Blakely_Figs_4_9_4_10';
            end

            switch type

                case 'Blakely_Figs_4_9_4_10'

                    % Replicate Blakely (1995) Figures 4.9 and 4.10 but
                    % using a small box rather than a dipole.
                    w = 1; V = w^3;
                    vx = [ -w/2 w/2 ];
                    vy = [ -w/2 w/2 ];
                    vz = [ -w/2 w/2 ];
                    vertBox = MagBox( vx, vy, vz, [ 0 0 1 ]'/V );
                    horizBox = MagBox( vx, vy, vz, [ 1 0 0 ]'/V );

                    % It will also be useful to define the dipole cases,
                    % just as a point of reference.
                    vertDipole = MagDipoles( [ 0 0 0 ]', [ 0 0 1 ]' );
                    horizDipole = MagDipoles( [ 0 0 0 ]', [ 1 0 0 ]' );

                    % Define a 2D plane for evaluation and show B field on that plane.
                    % This corresponds to Blakely Figure 4.9.
                    survey_plane = SurveyField( linspace(-2,2), linspace(-2,2), -1 );
                    fa = vertBox.showBfieldContours( 'z', survey_plane, 'z_down', true, 'view', [ -90 90 ], 'clim', [ -100 100 ], 'cbar', false, 'axlbl', '(a)' );
                    fc = vertBox.showBfieldContours( 'x', survey_plane, 'z_down', true, 'view', [ -90 90 ], 'clim', [ -100 100 ], 'cbar', false, 'axlbl', '(c)' );
                    fb = horizBox.showBfieldContours( 'x', survey_plane, 'z_down', true, 'view', [ -90 90 ], 'clim', [ -100 100 ], 'cbar', false, 'axlbl', '(b)' );
                    fd = horizBox.showBfieldContours( 'z', survey_plane, 'z_down', true, 'view', [ -90 90 ], 'clim', [ -100 100 ], 'cbar', false, 'axlbl', '(d)' );
                    fh49 = BaseTools.tileFigures( { fa fb fc fd }, 2, 2 );
                    
                    % Define a 1D transect and show B field along it.
                    % This corresponds to Blakely Figure 4.10.
                    survey_transect = SurveyField( linspace(-2,2), 0, -1 );
                    f_ca = vertBox.showBfieldContours( 'xz', survey_transect, 'ref_src', vertDipole, 'transect_labels', { '(c)' '(a)' } );
                    f_bd = horizBox.showBfieldContours( 'xz', survey_transect, 'ref_src', horizDipole, 'transect_labels', { '(b)' '(d)' } );
                    fh410 = BaseTools.tileFigures( { f_ca f_bd }, 1, 2 );
                    fh410.Position(3:4) = [ 1100 480 ];
                    
                    % Define a 3D region and show B field in that region.
                    survey_volume = SurveyField( linspace(-3,3), linspace(-3,3), linspace(-3,3) );
                    vfh = vertBox.showBfieldContours( 'xyz', survey_volume, 'z_down', true, 'view', [ -130 20 ], 'title', true );
                    hfh = horizBox.showBfieldContours( 'xyz', survey_volume, 'z_down', true, 'view', [ -130 20 ], 'title', true );

                    % For context, show the 1D and 2D surveys on top of the
                    % 3D volume plots.
                    survey_plane.overlaySurvey( vfh, 'clr', [ .5 0 .5 ] );
                    survey_plane.overlaySurvey( hfh, 'clr', [ .5 0 .5 ] );
                    survey_transect.overlaySurvey( vfh, 'clr', [ .9 0 .9 ] );
                    survey_transect.overlaySurvey( hfh, 'clr', [ .9 0 .9 ] );

                    % Assemble separate panels into a single figure.
                    fh3Dv = BaseTools.tileFigures( vfh, 1, 3 );
                    fh3Dh = BaseTools.tileFigures( hfh, 1, 3 );

                    % Set output arguments.
                    varargout{1} = { fh49, fh410, fh3Dv, fh3Dh };

                case 'HemingwayTikoo2018'

                    % Replicate figures from Hemingway & Tikoo (2018).
                    d = 0.5;
                    w = 1;
                    h = 2;
                    vx = [ -w/2 w/2 ];
                    vy = [ -Inf Inf ];
                    vz = [ -(d+h) -d ];
                    vertBox = MagBox( vx, vy, vz, [ 0 0 1 ]' );
                    horizBox = MagBox( vx, vy, vz, [ 1 0 0 ]' );

                    % Show field on 2D plane.
                    survey_plane = SurveyField( linspace(-3,3), 0, linspace(-3,1) );
                    vfh = vertBox.showBfieldContours( 'xz', survey_plane );
                    hfh = horizBox.showBfieldContours( 'xz', survey_plane );
                    BaseTools.tileFigures( [ hfh' vfh' ] );

                case 'askew_prism'

                    % Define a prism that is a tad askew.
                    roll = 10;
                    pitch = -20;
                    yaw = 30;
                    R_AS = BaseTools.rpy2rot( roll, pitch, yaw );
                    v_AS = [ 2 -.5 0 ]';
                    M = [ 0 1 0 ]'; % This is in the source coordinate frame.
                    myBox = MagBox( [ 0 2 ], [ 0 1 ], [ 0 0.6 ], M, v_AS, R_AS, 'S' );
                    myBox.drawBox( 'axes', true );

                    % Show field structure around this prism.
                    survey_plane = SurveyField( linspace(-8,8), linspace(-8,8), 0.5 );
                    BaseTools.tileFigures( myBox.showBfieldContours( 'xyz', survey_plane ) );
                    BaseTools.tileFigures( myBox.showBfieldContours( 'xyz' ) );
                    
                case 'singularities'

                    % Common settings.
                    M = [ 2 0.5 -1 ]'; % Keep the same magnetization for all.
                    d = 0.5; % Coarse grid spacing.
                    survey_volume = SurveyField( linspace(-10,10), linspace(-10,10), linspace(-10,10) );
                    z0 = 2; % Finite end position for semi-infinite prisms.

                    %
                    % Finite prism.
                    %
                    myBox = MagBox( [ 1 2 ], [ 2 4 ], [ 2 5 ], M );
                    BaseTools.tileFigures( myBox.showBfieldContours('xyz') );
                    xv = myBox.xrng(1)-1:d:myBox.xrng(2)+1;
                    yv = myBox.yrng(1)-1:d:myBox.yrng(2)+1;
                    zv = myBox.zrng(1)-1:d:myBox.zrng(2)+1;
                    myBox.showBfieldVectors( SurveyField( xv, yv, min(zv) ) );
                    myBox.showBfieldVectors( SurveyField( xv, min(yv), zv ) );
                    myBox.showBfieldVectors( SurveyField( min(xv), yv, zv ) );
                    myBox.showBfieldVectors( SurveyField( xv, yv, max(zv) ) );
                    myBox.showBfieldVectors( SurveyField( xv, max(yv), zv ) );
                    myBox.showBfieldVectors( SurveyField( max(xv), yv, zv ) );

                    %
                    % Semi-infinite (and quasi-semi-infinite) prisms.
                    %
                    for s = [ -1 1 ] % One side and then the other.
                        for far = [ 100 Inf ] % For each side, do a finite and infinite limit for the far end.
                            myBox = MagBox( [ 1 2 ], [ 2 4 ], sign(s)*[ z0 far ], M );
                            BaseTools.tileFigures( myBox.showBfieldContours('xyz',survey_volume) );
                            myBox.showBfieldVectors( SurveyField( xv, yv, 0 ) );
                        end
                    end

                    %
                    % Infinite prism.
                    %
                    myBox = MagBox( [ 1 2 ], [ 2 4 ], [ -Inf Inf ], M );
                    BaseTools.tileFigures( myBox.showBfieldContours('xyz') );

                    %
                    % Wide slab.
                    %
                    wide_survey = SurveyField( linspace(-300,300), linspace(-300,300), linspace(-200,200) );
                    myBox = MagBox( [ -100 100 ], [ -200 200 ], [ 2 5 ], M );
                    BaseTools.tileFigures( myBox.showQfieldContours('all',wide_survey) );
                    myBox = MagBox( [ -Inf Inf ], [ -200 200 ], [ 2 5 ], M );
                    BaseTools.tileFigures( myBox.showQfieldContours('all',wide_survey) );
                    myBox = MagBox( [ -Inf Inf ], [ -Inf 200 ], [ 2 5 ], M );
                    BaseTools.tileFigures( myBox.showQfieldContours('all',wide_survey) );                   
                    myBox = MagBox( [ -Inf 100 ], [ -Inf 200 ], [ 2 5 ], M );
                    BaseTools.tileFigures( myBox.showQfieldContours('all',wide_survey) );                   
                    myBox = MagBox( [ -Inf 100 ], [ -Inf 200 ], [ -Inf 5 ], M );
                    BaseTools.tileFigures( myBox.showQfieldContours('all',wide_survey) );                   

                case 'staticQ'

                    % Construct a prism with static dimensions but then
                    % shift the magnetization direction. Here, we're
                    % testing to see how performance is improved if we
                    % avoid recalculating Q every time.
                    myBox = MagBox( [ -1 3 ], [ 1 5 ], [ 1 2 ], [ 1 0 0 ]' );

                    % Define a survey transect.
                    survey = SurveyField( linspace(-4,4), 0, 0 );

                    % Let's have M vary sinusoidally.
                    theta = linspace(0,2*pi,361);
                    M = [ cos(theta); sin(theta); zeros(size(theta)) ];
                    num_thetas = length(theta);

                    % Now define a B grid that runs across the survey transect
                    % and has as many rows as there are M orientations.
                    % Let's consider only Bx for simplicity.
                    Bx_grid = nan(num_thetas,survey.Np);

                    % Populate the B grid.
                    tic;
                    for k = 1 : num_thetas
                        myBox.Mr = M(:,k);
                        thisB = myBox.computeBfield(survey.p);
                        Bx_grid(k,:) = thisB(1,:);
                    end
                    t_elapsed = toc;
                    fprintf( 'Computation time without pre-computing Q = %s\n', BaseTools.timetostr(t_elapsed) );

                    % Show results.
                    figure;
                    imagesc( survey.xv, theta, Bx_grid );
                    colormap( MagSource.magcolors );
                    colorbar;
                    xlabel( 'x' );
                    ylabel( '\theta' );

                    % Do the same thing again, but this time compute Q only
                    % once.
                    Q = myBox.computeQfield(survey.p);
                    Bx_grid = nan(num_thetas,survey.Np);
                    tic;
                    for k = 1 : num_thetas
                        myBox.Mr = M(:,k);
                        thisB = myBox.computeBfield(survey.p,'Q',Q);
                        Bx_grid(k,:) = thisB(1,:);
                    end
                    t_elapsed = toc;
                    fprintf( 'Computation time with pre-computing Q = %s\n', BaseTools.timetostr(t_elapsed) );

                    % Show results.
                    figure;
                    imagesc( survey.xv, theta, Bx_grid );
                    colormap( MagSource.magcolors );
                    colorbar;
                    xlabel( 'x' );
                    ylabel( '\theta' );

            end            

        end

    end
    
end




