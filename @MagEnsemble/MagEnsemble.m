classdef MagEnsemble < MagSource
    
    % MagEnsemble: class representing an array of magnetic source bodies.
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
        
        sources % array of magnetic sources of possibly different types (e.g., MagBox, MagDipoles).
        
    end
    
    properties (Dependent)

        mtot % total magnetic moment (Am^2)
        Vtot % total volume represented by the array of sources
        M % magnetic moment per unit volume (A/m), may or may not vary between dipoles.
        xrng % (1x2 array) range of x positions spanned by the sources.
        yrng % (1x2 array) range of y positions spanned by the sources.
        zrng % (1x2 array) range of z positions spanned by the sources.
        N % number of sources in the array.
        
    end
    
    methods
        
        % Constructor.
        function newSource = MagEnsemble( new_sources )
            newSource.sources = new_sources;
        end

        % Simple getter methods.
        function mtot = get.mtot( thisMagEnsemble )
            mtot = 0;
            for i = 1 : length(thisMagEnsemble.sources)
                mtot = mtot + thisMagEnsemble.sources(i).mtot;
            end
        end
        function Vtot = get.Vtot( thisMagEnsemble )
            Vtot = 0;
            for i = 1 : length(thisMagEnsemble.sources)
                Vtot = Vtot + thisMagEnsemble.sources(i).Vtot;
            end
        end
        function M = get.M( thisMagEnsemble )
            M = thisMagEnsemble.mtot / thisMagEnsemble.Vtot;
        end
        function xrng = get.xrng( thisMagEnsemble )
            xvals = nan(length(thisMagEnsemble.sources),2);
            for i = 1 : length(thisMagEnsemble.sources)
                xvals(i,:) = thisMagEnsemble.sources(i).xrng;
            end
            xrng = [ min(xvals(:)) max(xvals(:)) ];
        end
        function yrng = get.yrng( thisMagEnsemble )
            yvals = nan(length(thisMagEnsemble.sources),2);
            for i = 1 : length(thisMagEnsemble.sources)
                yvals(i,:) = thisMagEnsemble.sources(i).yrng;
            end
            yrng = [ min(yvals(:)) max(yvals(:)) ];
        end
        function zrng = get.zrng( thisMagEnsemble )
            zvals = nan(length(thisMagEnsemble.sources),2);
            for i = 1 : length(thisMagEnsemble.sources)
                zvals(i,:) = thisMagEnsemble.sources(i).zrng;
            end
            zrng = [ min(zvals(:)) max(zvals(:)) ];
        end
        function N = get.N( thisMagEnsemble )
            N = length(thisMagEnsemble.sources);
        end
        
        % Compute fields.
        function B = computeBfield( thisMagEnsemble, p, varargin )
            Np = size(p,2);
            B = zeros(3,Np);
            for k = 1 : thisMagEnsemble.N
                B = B + thisMagEnsemble.sources(k).computeBfield( p, varargin{:} );
            end
        end
        function Q = computeQfield( thisMagEnsemble, p, varargin )
            Np = size(p,2);
            Q = zeros(3,3,Np);
            for k = 1 : thisMagEnsemble.N
                Q = Q + thisMagEnsemble.sources(k).computeQfield( p, varargin{:} );
            end
        end

        % Display functions.
        function drawEnsemble( thisMagEnsemble, varargin )
            if ~isempty(varargin) && ishandle(varargin{1})
                ah = varargin{1};
                opt_args = varargin(2:end);
            else
                fh = figure;
                ah = axes('Parent',fh); 
                xlabel('x');
                ylabel('y');
                zlabel('z');
                opt_args = varargin;
            end
            hold(ah,'on');
            for k = 1 : thisMagEnsemble.N
                thisMagEnsemble.sources(k).drawSource( ah, opt_args{:} );
            end
        end
        
    end

    methods (Static)

        % Testing.
        function varargout = unit_test( type )

            close all;
            if ~exist( 'type', 'var' )
                type = 'baseline';
            end
            switch type

                case 'baseline'

                    % Create source array.
                    myEnsemble = MagEnsemble( [ ...
                        MagBox( [ -4.5 0 ], [ -1 5 ], [ -0.5 -1 ], [ 1 2 0 ]' ); ...
                        MagBox( [ -3 -2.5 ], [ -5 1 ], [ -2 -5 ], [ 10 0 0 ]' ); ...
                        MagBox( [ 3 4 ], [ -3 4 ], [ -1 -3 ], [ 5 5 0 ]' ); ...
                        MagDipoles( linspace( 1, 2, 3 ), linspace( -1, 2, 5 ), linspace( -1, -2, 5 ), [ 0 5 1 ]' ) ...
                        ] );

                    % Show field on a couple of 2D planes.
                    survey_planes = [
                        SurveyField( linspace(-5,5), linspace(-5,5), 1 )
                        SurveyField( linspace(-5,5), 1, linspace(-6,1) )
                        ];                    
                    BaseTools.tile_figures( myEnsemble.showBfieldContours( 'xyz', survey_planes, 'view', [ 25 15 ] ) );
                    BaseTools.tile_figures( myEnsemble.showQfieldContours( { 'xx' 'xy' }, survey_planes, 'view', [ 25 15 ] ) );

                    % Show field over 3D volume.
                    survey_volume = SurveyField( linspace(-10,10), linspace(-10,10), linspace(0,10) );
                    myEnsemble.showBfieldContours( 'z', survey_volume );
                    BaseTools.tile_figures( myEnsemble.showQfieldContours( 'all', survey_volume, 'view', [ 25 15 ] ) );

                case 'Bongiolo4'

                    % Need to translate from Bongiolo's Table 1 into my
                    % terminology and then compare to their Figure 4.
                    inc = 5; % degrees
                    dec = 10; % degrees
                    Bint = 439.82e-9; % Teslas
                    B_ambient = MagSource.B_inc_dec_to_enr( Bint, inc, dec );
                    X = 1;
                    Mr = 0;
                    wx = 20;
                    wy = 20;
                    wz = 2;
                    x0 = 30;
                    y0 = 30;
                    topdepth = 1;
                    yinc = 0;
                    myBox = BongioloMagBox( B_ambient, X, Mr, wx, wy, wz, x0, y0, topdepth, yinc );
                    myBox.Bref = B_ambient;
                    survey_plane = SurveyField( linspace(0,64), linspace(0,64), 0 );
                    myBox.showBfieldContours( 'a', survey_plane, 'view', [ 0 90 ] );
                    clim( [ -44 39 ] ); colormap(jet);
                    BaseTools.tile_figures( myBox.showBfieldContours( 'xyz', survey_plane, 'view', [ 0 90 ] ) );
                    survey_volume = SurveyField( linspace(0,64), linspace(0,64), linspace(-20,20) );
                    BaseTools.tile_figures( myBox.showBfieldContours( 'xyz', survey_volume, 'view', [ 0 90 ] ) );

                case 'Bongiolo7'

                    % Need to translate from Bongiolo's Table 3 into my
                    % terminology and then compare to their Figure 7.
                    dec = -13.3; % degrees
                    inc = 12.56667; % degrees
                    Bint = 27865e-9; % Teslas
                    B_ambient = MagSource.B_inc_dec_to_enr( Bint, inc, dec );
                    myEnsemble = MagEnsemble( [ ...
                        BongioloMagBox( B_ambient, 0.027, 0, 2000, 200, 500, 3500, 3500, 80, 25 ) ...
                        BongioloMagBox( B_ambient, 0.027, 0, 2000, 200, 500, 1500, 1500, 50, -25 ) ...
                        BongioloMagBox( B_ambient, 0.027, 0, 500, 500, 500, 4500, 4500, 200, 0 ) ...
                        BongioloMagBox( B_ambient, 0.027, 0, 1000, 500, 500, 1500, 4500, 100, 45 ) ...
                        BongioloMagBox( B_ambient, 0.027, 0, 1500, 100, 250, 4500, 1000, 50, -215 ) ...
                        BongioloMagBox( B_ambient, 0.027, 0, 1000, 200, 500, 4500, 1000, 150, 75 ) ...
                        ] );
                    myEnsemble.Bref = B_ambient;
                    survey_plane = SurveyField( linspace(0,6000), linspace(0,6000), 0 );
                    myEnsemble.showBfieldContours( 'a', survey_plane, 'view', [ 0 90 ] );
                    clim( [ -275 95 ] ); colormap(jet);
                    BaseTools.tile_figures( myEnsemble.showBfieldContours( 'xyz', survey_plane, 'view', [ 0 90 ] ) );

            end

        end

    end
     
end




function thisBox = BongioloMagBox( B_ambient, X, Mr, wx, wy, wz, x0, y0, topdepth, yinc )

% Translate from Bongiolo+ terminology to mine.

% Note: Bongiolo+ seems to be using a left-handed coordinate system with +y
% aligned with geographic north, +x aligned with geographic east, and +z
% pointing downward. I always use a right-handed system, so I will adopt 
% +x (east), +y (north), +z (upward).

% Compute induced rock magnetization given the specified magnetic
% susceptibility and the ambient B-field.
M_ind = B_ambient * X/MagSource.mu_naught;

% If a remanent magnetization has also been specified, add that here too.
M = M_ind + Mr;

% Bongiolo+ allow for a prism to be centered away from the origin, and
% rotated in the xy plane (rotation about the z-axis). We have to be
% careful about the sign of the rotation here. 
v_AS = [ x0 y0 -(topdepth+wz/2) ]';
R_AS = BaseTools.rpy2rot( 0, 0, yinc )';

% Finally, we have the necessary information to specify the prism
% dimensions and effective magnetization. The last argument tells MagBox
% that the magnetization direction is specified in the analysis ('A')
% coordinate system and not a coordinate system that is rotated with the
% prism.
thisBox = MagBox( [ -wx/2 wx/2 ], [ -wy/2 wy/2 ], [ -wz/2 wz/2 ], M, v_AS, R_AS, 'A' );

end


