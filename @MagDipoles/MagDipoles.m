classdef MagDipoles < MagSource
    
    % MagDipoles: class representing a magnetic dipole or an array of
    % dipoles.
    %
    %   MagDipoles( q, m ) is the standard constructor call, where
    %       q is a 3xN array representing the N dipole locations and
    %       m is a 3xN array representing the magnetic moment of each dipole
    %
    %   MagDipoles( xv, yv, zv, m ) is an alternate form, where
    %       xv, yv, zv are arrays of x, y, z positions used to generate a
    %       meshgrid of dipoles filling the space defined by these points.
    %       m is a 3x1 vector representing the total magnetic moment (this
    %       will be divided evenly among the array of dipoles). The spacing
    %       between the dipoles will be used to compute the "effective
    %       volume" represented by each dipole.
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
    
    properties (Access=public)
        
        q % (3xN vector) position of the N dipoles (meters)
        m % (3xN vector) magnetic moment of each dipole (Am^2)
        dV % (1xN vector) volume corresponding to each dipole (m^3)
        
    end
    
    properties (Dependent)
        
        mtot % total magnetic moment (Am^2)
        Vtot % total volume represented by the array of dipoles (N*dV)
        M % magnetic moment per unit volume (A/m), may or may not vary between dipoles.
        Nq % number of dipoles in the array
        xrng % (1x2 array) range of x positions spanned by the dipoles.
        yrng % (1x2 array) range of y positions spanned by the dipoles.
        zrng % (1x2 array) range of z positions spanned by the dipoles.
        
    end
    
    methods
        
        % Constructor.
        function newDipoles = MagDipoles( varargin )
            % TO DO: make a better constructor.           
            if nargin == 1 && ischar( varargin{1} ) && strcmp( varargin{1}, 'demo' )
                % Here, the user wants a demo.
                newDipoles.q = [ 3 4 0 ]';
                newDipoles.m = [ 1 0 0 ]';
            elseif nargin == 2 || nargin == 3
                % Here, we have the q and m arrays specified directly.
                newDipoles.q = varargin{1};
                newDipoles.m = varargin{2};
                if nargin == 3
                    newDipoles.dV = varargin{3};
                else
                    newDipoles.dV = 1;
                end
            elseif nargin == 4 || nargin == 5
                % Here, the user is supplying the dimensions of a box,
                % which we will fill with dipoles at the specified
                % resolution and with the specified equivalent total
                % magnetic moment (specified as a 3x1 vector, M).
                xvals = varargin{1};
                yvals = varargin{2};
                zvals = varargin{3};
                m = varargin{4};
                [ qx, qy, qz ] = meshgrid( xvals, yvals, zvals );
                N = numel(qx);
                newDipoles.q = [ qx(:)'; qy(:)'; qz(:)' ];
                newDipoles.m = repmat( m/N, 1, N );
                if nargin == 5
                    newDipoles.dV = varargin{5};
                else
                    newDipoles.dV = 1;
                end
            end
        end

        % Simple getter methods.
        function mtot = get.mtot( theseDipoles )
            mtot = sum( sqrt( theseDipoles.m(1,:).^2 + theseDipoles.m(2,:).^2 + theseDipoles.m(3,:).^2 ), 2 );
        end
        function Vtot = get.Vtot( theseDipoles )
            if isscalar(theseDipoles.dV)
                Vtot = theseDipoles.Nq * theseDipoles.dV;
            else
                Vtot = sum(theseDipoles.dV);
            end
        end
        function M = get.M( theseDipoles )
            if isempty(theseDipoles.dV)
                error( 'dV must be specified before we can compute M.' );
            else
                M = theseDipoles.m ./ theseDipoles.dV;
                if all(all(M == M(:,1)))
                    M = M(:,1);
                end
            end
        end
        function Nq = get.Nq( theseDipoles )
            [ qrows, qcols ] = size( theseDipoles.q );
            if qrows ~= 3
                error( 'q must have 3 rows' );
            end
            Nq = qcols;
        end
        function xrng = get.xrng( theseDipoles )
            xrng = [ min(theseDipoles.q(1,:)) max(theseDipoles.q(1,:)) ];
        end
        function yrng = get.yrng( theseDipoles )
            yrng = [ min(theseDipoles.q(2,:)) max(theseDipoles.q(2,:)) ];
        end
        function zrng = get.zrng( theseDipoles )
            zrng = [ min(theseDipoles.q(3,:)) max(theseDipoles.q(3,:)) ];
        end

        % Setter methods.
        function theseDipoles = set.M( theseDipoles, new_M )
            if isempty(theseDipoles.dV)
                error( 'dV must be specified before we can set M.' );
            else
                % Here, each dipole has an associated volume, dV. Use that
                % to compute the magnetic moment we need to assign to each.
                theseDipoles.m = new_M * theseDipoles.dV;
            end
        end
        
        % Compute fields.
        Q = computeQfield( theseDipoles, p, varargin )
        B = computeBfield( theseDipoles, p, varargin )

        % Display functions.
        function ah = drawDipoles( theseDipoles, varargin )
            if ~isempty(varargin) && ishandle(varargin{1})
                ah = varargin{1};
                opt_args = varargin(2:end);
            else
                opt_args = varargin;
            end
            def_args = { 'overlay', false, 'axislabels', 'xyz', 'clr', [ .2 .2 .2 ] };
            args = BaseTools.argarray2struct( [ def_args opt_args{:} ] );
            if ~exist( 'ah', 'var' ) || ~ishandle( ah )
                fh = figure;
                ah = axes('Parent',fh);
                xlabel('x');
                ylabel('y');
                zlabel('z');
            end
            hold(ah,'on');
            grid(ah,'on');
            axis(ah,'tight');
            corners = SurveyField( theseDipoles.xrng, theseDipoles.yrng, theseDipoles.zrng );
            if sum(corners.nonsingdims)==3
                axis(ah,'equal'); % This array of dipoles spans three dimensions so all the plot's axes will be spatial.
            end
            for k = 1 : theseDipoles.Nq
                plot3( theseDipoles.q(1,k), theseDipoles.q(2,k), theseDipoles.q(3,k), 'o', 'Color', args.clr, 'MarkerSize', 9, 'LineWidth', 2.0 );
            end
        end

    end
    
    methods (Static)

        function varargout = unit_test( type )

            close all;
            if ~exist( 'type', 'var' )
                type = 'Blakely_Figs_4_9_4_10';
            end

            switch type

                case 'Blakely_Figs_4_9_4_10'

                    % Define two dipoles, one magnetized vertically and one
                    % horizontally.
                    vertDipole = MagDipoles( [ 0 0 0 ]', [ 0 0 1 ]' );
                    horizDipole = MagDipoles( [ 0 0 0 ]', [ 1 0 0 ]' );

                    % Define a 2D plane for evaluation and show B field on that plane.
                    % This corresponds to Blakely Figure 4.9.
                    survey_plane = SurveyField( linspace(-2,2), linspace(-2,2), -1 );
                    fa = vertDipole.showBfieldContours( 'z', survey_plane, 'view', [ 90 -90 ], 'clim', [ -100 100 ], 'cbar', false, 'axlbl', '(a)' );
                    fc = vertDipole.showBfieldContours( 'x', survey_plane, 'view', [ 90 -90 ], 'clim', [ -100 100 ], 'cbar', false, 'axlbl', '(c)' );
                    fb = horizDipole.showBfieldContours( 'x', survey_plane, 'view', [ 90 -90 ], 'clim', [ -100 100 ], 'cbar', false, 'axlbl', '(b)' );
                    fd = horizDipole.showBfieldContours( 'z', survey_plane, 'view', [ 90 -90 ], 'clim', [ -100 100 ], 'cbar', false, 'axlbl', '(d)' );
                    fh49 = BaseTools.tile_figures( { fa fb fc fd }, 2, 2 );
                    
                    % Define a 1D transect and show B field along it.
                    % This corresponds to Blakely Figure 4.10.
                    survey_transect = SurveyField( linspace(-2,2), 0, -1 );
                    f_ca = vertDipole.showBfieldContours( 'xz', survey_transect, 'transect_labels', { '(c)' '(a)' } );
                    f_bd = horizDipole.showBfieldContours( 'xz', survey_transect, 'transect_labels', { '(b)' '(d)' } );
                    fh410 = BaseTools.tile_figures( { f_ca f_bd }, 1, 2 );
                    fh410.Position(3:4) = [ 1100 480 ];
                    
                    % Define a 3D region and show B field in that region.
                    survey_volume = SurveyField( linspace(-3,3), linspace(-3,3), linspace(-3,3) );
                    vfh = vertDipole.showBfieldContours( 'xyz', survey_volume, 'z_down', true, 'view', [ -130 -20 ] );
                    hfh = horizDipole.showBfieldContours( 'xyz', survey_volume, 'z_down', true, 'view', [ -130 -20 ] );

                    % For context, show the 1D and 2D surveys on top of the
                    % 3D volume plots.
                    survey_plane.overlaySurvey( vfh, 'clr', [ .5 0 .5 ] );
                    survey_plane.overlaySurvey( hfh, 'clr', [ .5 0 .5 ] );
                    survey_transect.overlaySurvey( vfh, 'clr', [ .9 0 .9 ] );
                    survey_transect.overlaySurvey( hfh, 'clr', [ .9 0 .9 ] );

                    % Assemble separate panels into a single figure.
                    fh3Dv = BaseTools.tile_figures( vfh, 1, 3 );
                    fh3Dh = BaseTools.tile_figures( hfh, 1, 3 );

                    % Set output arguments.
                    varargout{1} = { fh49, fh410, fh3Dv, fh3Dh };

                    % Just out of interest, let's also look at Q.
                    % Recall that Q is fundamentally related to the volume
                    % occupied by the source, and so some volume has to be
                    % assumed here. By default, the assumption is that each
                    % dipole corresponds to a volume of 1x1x1 length units.
                    vertDipole.showQfieldContours( 'all', survey_transect );
                    BaseTools.tile_figures( vertDipole.showQfieldContours( 'all', survey_plane ), 3, 3 );
                    BaseTools.tile_figures( vertDipole.showQfieldContours( 'all', survey_volume, 'z_down', true, 'view', [ -130 -20 ] ), 3, 3 );
                    
                case 'multiple_dipoles'

                    % Define array of dipoles.
                    q = [ 
                        0 3
                        0 2
                        0 0
                        ];
                    m = [ 3 .1 -1 ]';
                    myDipoles = MagDipoles( q, m );
                    b = 6;
                                        
                    % Define a plane for evaluation and show Q field on that plane.
                    survey_plane = SurveyField( linspace(-b,b), linspace(-b,b), -0.1 );
                    BaseTools.tile_figures( myDipoles.showQfieldContours( 'all', survey_plane ), 3, 3 );
                    BaseTools.tile_figures( myDipoles.showBfieldContours( 'xyz', survey_plane ), 1, 3 );

                    % Define a survey volume and show field contours.
                    survey_volume = SurveyField( linspace(-b,b), linspace(-b,b), linspace(-b,b) );
                    BaseTools.tile_figures( myDipoles.showQfieldContours( 'all', survey_volume ), 3, 3 );
                    BaseTools.tile_figures( myDipoles.showBfieldContours( 'xyz', survey_volume ), 1, 3 );

                case 'dipole_grid'

                    % Define grid.
                    xv = linspace( 1, 3, 5 );
                    yv = linspace( 2, 5, 4 );
                    zv = linspace( 1, 1.5, 3 );
                    myDipoles = MagDipoles( xv, yv, zv, [ 3 .1 -1 ]' );
                    myDipoles.drawDipoles;
                    b = 12;

                    % Define a plane for evaluation and show Q field on that plane.
                    survey_plane = SurveyField( linspace(-b,b), linspace(-b,b), -1 );
                    BaseTools.tile_figures( myDipoles.showQfieldContours( 'all', survey_plane ), 3, 3 );
                    BaseTools.tile_figures( myDipoles.showBfieldContours( 'xyz', survey_plane ), 1, 3 );
                    
                    % Define a survey volume and show field contours.
                    n = 50; % For 3D volumes, it may be useful to limit resolution.
                    survey_volume = SurveyField( linspace(-b,b,n), linspace(-b,b,n), linspace(-b,b,n) );
                    BaseTools.tile_figures( myDipoles.showQfieldContours( 'all', survey_volume ), 3, 3 );
                    BaseTools.tile_figures( myDipoles.showBfieldContours( 'xyz', survey_volume ), 1, 3 );
                    
                    % % Show field contours.
                    % BaseTools.tile_figures( myDipoles.showQfieldContours( 'all', [] ), 3, 3 );
                    % BaseTools.tile_figures( myDipoles.showBfieldContours( 'xyz', [] ), 1, 3 );

                case 'timing'

                    % We will build grids of dipoles and see how long it
                    % takes to compute the B-field. We'll vary both the
                    % density of the grid and the density of the survey.
                    Nsmooth = 20;
                    grid_res = 3:6;
                    survey_res = [ 20 30 40 ];
                    M = [ 3 .1 -1 ]';
                    Nq = grid_res.^3;
                    Np = survey_res.^3;
                    Nevals = Nq' * Np;

                    % Establish an array of output times.
                    t = zeros( length(grid_res), length(survey_res) );

                    % Consider each combination in turn.
                    fprintf( 'Starting timing test. Please stand by...\n' );
                    for i = 1 : length(grid_res)

                        % Build dipole array at this resolution.
                        myDipoles = MagDipoles( linspace(-1,1,grid_res(i)), linspace(-1,1,grid_res(i)), linspace(-1,1,grid_res(i)), M );

                        % Check how long it takes for each of several
                        % different survey resolutions.
                        for j = 1 : length(survey_res)

                            % Build a 3D survey at this resolution.
                            survey = SurveyField( linspace(-10,10,survey_res(j)), linspace(-10,10,survey_res(j)), linspace(-10,10,survey_res(j)) );

                            % Compute B-field (take the average of Nsmooth calculations).
                            tic;
                            for k = 1 : Nsmooth
                                myDipoles.computeBfield(survey.p);
                            end
                            t(i,j) = toc/Nsmooth;

                        end

                    end

                    % Show results graphically.
                    figure;
                    imagesc( Np, Nq, t );
                    xlabel( 'Np' );
                    ylabel( 'Nq' );
                    colorbar;

                    % Now normalize by Nevals.
                    tau = t./Nevals;
                    figure;
                    imagesc( Np, Nq, tau );
                    xlabel( 'Np' );
                    ylabel( 'Nq' );
                    colorbar;
                    tau_mean = mean( tau(:) );

                    % Show tau vs Nevals.
                    figure;
                    plot( Nevals, t, 'LineWidth', 2.0 );
                    xlabel( 'Np * Nq' );
                    ylabel( 'Computation Time (s)' );
                    grid on;
                    fprintf( 'Mean computation time: %f ns per evaluation.\n', tau_mean*1e9 );                    
                    fprintf( 'Evaluations per second: %e.\n', 1/tau_mean );                    

                case 'field_lines'

                    % Define a dipole.
                    myDipole = MagDipoles( [ 0 0 0 ]', [ 0 0 1 ]' );

                    % Define a survey region.
                    survey_volume = SurveyField( linspace(-3,3), linspace(-3,3), linspace(-3,3) );

                    % Show field lines in this region.
                    myDipole.showBfieldLines( survey_volume );

            end

        end

    end

end



