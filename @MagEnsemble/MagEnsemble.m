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
        function ah = drawEnsemble( thisMagEnsemble, varargin )
            [ args, opt_args ] = BaseTools.argarray2struct( varargin );
            ah = BaseTools.extractAxesHandle( args );
            xlabel('x');
            ylabel('y');
            zlabel('z');
            hold(ah,'on');
            for k = 1 : thisMagEnsemble.N
                thisMagEnsemble.sources(k).drawSource( ah, opt_args{:} );
            end
        end
        
    end

    methods (Static)

        % Testing.
        function varargout = unit_test( type )

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
                    BaseTools.tileFigures( myEnsemble.showBfieldContours( 'xyz', survey_planes, 'view', [ 25 15 ] ) );
                    BaseTools.tileFigures( myEnsemble.showQfieldContours( { 'xx' 'xy' }, survey_planes, 'view', [ 25 15 ] ) );

                    % Show field over 3D volume.
                    survey_volume = SurveyField( linspace(-10,10), linspace(-10,10), linspace(0,10) );
                    myEnsemble.showBfieldContours( 'z', survey_volume );
                    BaseTools.tileFigures( myEnsemble.showQfieldContours( 'all', survey_volume, 'view', [ 25 15 ] ) );

                case 'mixed_orientation'

                    % Create source array.
                    myEnsemble = MagEnsemble( [ ...
                        MagBox( [ -4 0 ], [ -1 3 ], [ -1 1 ], [ 1 0 0 ]' ); ...
                        MagBox( [ 0 5 ], [ 0 2 ], [ 0 1 ], [ 2 0 0 ]', [ 4 0 -2 ]', BaseTools.rpy2rot( 17, -28, 39 ), 'S' ); ...
                        ] );

                    % Show field structure around this prism.
                    survey_plane = SurveyField( linspace(-8,8), linspace(-8,8), 0 );
                    fh = myEnsemble.showBfieldContours( 'x', survey_plane, 'view', [ 20 40 ] );
                    fh.CurrentAxes.ZLim = [ -2 2 ];
                    fh.Position(4) = 400;
                    varargout{1} = fh;

                    % Show composite Q.
                    BaseTools.tileFigures( myEnsemble.showQfieldContours( 'all' ) );
                                        
                case 'Bongiolo4'

                    % Need to translate from Bongiolo's Table 1 into my
                    % terminology and then compare to their Figure 4.
                    inc = 5; % degrees
                    dec = 10; % degrees
                    Bint = 439.82e-9; % Teslas
                    B_ambient = MagSource.B_inc_dec_to_enr( Bint, inc, dec );
                    chi = 1;
                    Mr = zeros(3,1);
                    wx = 20;
                    wy = 20;
                    wz = 2;
                    x0 = 30;
                    y0 = 30;
                    topdepth = 1;
                    yinc = 0;
                    myBox = BongioloMagBox( chi, Mr, wx, wy, wz, x0, y0, topdepth, yinc );
                    myBox.Bref = B_ambient;
                    survey_plane = SurveyField( linspace(0,64), linspace(0,64), 0 );
                    myBox.showBfieldContours( 'a', survey_plane, 'view', [ 0 90 ] );
                    BaseTools.tileFigures( myBox.showBfieldContours( 'xyz', survey_plane, 'view', [ 0 90 ] ) );
                    survey_volume = SurveyField( linspace(0,64), linspace(0,64), linspace(-20,20) );
                    BaseTools.tileFigures( myBox.showBfieldContours( 'xyz', survey_volume, 'view', [ 0 90 ] ) );

                case 'Bongiolo7'

                    % Need to translate from Bongiolo's Table 3 into my
                    % terminology and then compare to their Figure 7.
                    dec = -13.3; % degrees
                    inc = 12.56667; % degrees
                    Bint = 27865e-9; % Teslas
                    B_ambient = MagSource.B_inc_dec_to_enr( Bint, inc, dec );
                    myEnsemble = MagEnsemble( [ ...
                        BongioloMagBox( 0.027, zeros(3,1), 2000, 200, 500, 3500, 3500, 80, 25 ) ...
                        BongioloMagBox( 0.027, zeros(3,1), 2000, 200, 500, 1500, 1500, 50, -25 ) ...
                        BongioloMagBox( 0.027, zeros(3,1), 500, 500, 500, 4500, 4500, 200, 0 ) ...
                        BongioloMagBox( 0.027, zeros(3,1), 1000, 500, 500, 1500, 4500, 100, 45 ) ...
                        BongioloMagBox( 0.027, zeros(3,1), 1500, 100, 250, 4500, 1000, 50, -215 ) ...
                        BongioloMagBox( 0.027, zeros(3,1), 1000, 200, 500, 4500, 1000, 150, 75 ) ...
                        ] );
                    myEnsemble.Bref = B_ambient;
                    survey_plane = SurveyField( linspace(0,6000), linspace(0,6000), 0 );
                    fh_total_anom = myEnsemble.showBfieldContours( 'a', survey_plane, 'view', [ 0 90 ] );
                    fh_B_components = BaseTools.tileFigures( myEnsemble.showBfieldContours( 'xyz', survey_plane, 'view', [ 0 90 ] ) );

                    % Set output arguments.
                    varargout{1} = { fh_total_anom, fh_B_components };

            end

        end

    end
     
end




function thisBox = BongioloMagBox( user_chi, Mr, wx, wy, wz, x0, y0, topdepth, yinc )

% Translate from Bongiolo+ terminology to mine.

% Note: Bongiolo+ seems to be using a left-handed coordinate system with +y
% aligned with geographic north, +x aligned with geographic east, and +z
% pointing downward. 
% I will use a right-handed system with +x (east), +y (north), +z (upward). 

% Bongiolo+ allow for a prism to be centered away from the origin, and
% rotated in the xy plane (rotation about the z-axis). We have to be
% careful about the sign of the rotation here. 
v_AS = [ x0 y0 -(topdepth+wz/2) ]';
R_AS = BaseTools.rpy2rot( 0, 0, yinc )';

% Finally, we have the necessary information to specify the prism
% dimensions, remanent magnetization, and susceptibility. The total
% magnetization (remanent plus induced) will be calculated later when the B
% field is requested.
thisBox = MagBox( [ -wx/2 wx/2 ], [ -wy/2 wy/2 ], [ -wz/2 wz/2 ], Mr, v_AS, R_AS, Mframe='A', chi=user_chi );

end


