classdef SurveyField

    % SurveyField: simple class for holding an array of points constituting
    % some survey field. 
    %
    %   Standard constructor is: SurveyField( xv, yv, zv ), where xv, yv, zv
    %   are the ranges of x, y, z values that define a meshgrid.
    %
    %   Alternate: SurveyField( p ), where p is a 3xN array, allows for
    %   non-meshgrid layout of evaluation points.
    %
    %   Bonus: SurveyField( source ), where source is a MagSource object,
    %   automatically defines a survey field based on the size of the
    %   source.
    %
    %
    %   Disclaimer: This code is provided as-is, has been tested only very
    %   informally, and may not always behave as intended. I find it useful for
    %   my own work, and I hope you will too, but I make no guarantees as to
    %   the accuracy or robustness of this code. This code is also actively
    %   under development and future versions may not be backward compatible.
    %
    %   Created: 2025-06-07
    %   Doug Hemingway (douglas.hemingway@utexas.edu)
    %   University of Texas at Austin
    %
    
    properties (Access=private)
        int_xv % range of x values for a meshgrid
        int_yv % range of y values for a meshgrid
        int_zv % range of z values for a meshgrid
        int_X % x coordinates in a meshgrid
        int_Y % y coordinates in a meshgrid
        int_Z % z coordinates in a meshgrid
        int_p % 3xN array of N 3x1 evaluation point vectors
    end

    properties (Dependent)
        xv % range of x values for a meshgrid
        yv % range of y values for a meshgrid
        zv % range of z values for a meshgrid
        X % x coordinates in a meshgrid
        Y % y coordinates in a meshgrid
        Z % z coordinates in a meshgrid
        p % 3xN array of N 3x1 evaluation point vectors
    end

    properties (Dependent)
        Np % Number of points in the array
        dimsizes % Dimensions of the grid (if it is a grid).
        singdims % How many singleton dimensions are there?
        nonsingdims % How many non-singleton dimensions are there?
    end

    methods

        % Constructor.
        function survey = SurveyField( varargin )
            if nargin == 3
                % Here, the user has supplied the meshgrid coordinates.
                survey.int_xv = varargin{1};
                survey.int_yv = varargin{2};
                survey.int_zv = varargin{3};
                [ survey.int_X, survey.int_Y, survey.int_Z ] = meshgrid( survey.xv, survey.yv, survey.zv );
                survey.int_p = [ survey.X(:)'; survey.Y(:)'; survey.Z(:)' ];
            elseif nargin == 1
                if isa( varargin{1}, 'MagSource' )
                    % Here, we have a MagSource object. Define a survey
                    % that covers the source body plus some margin that
                    % depends on the scale of the body.
                    src = varargin{1};
                    wx = src.xrng(2)-src.xrng(1);
                    wy = src.yrng(2)-src.yrng(1);
                    wz = src.zrng(2)-src.zrng(1);
                    ws = [ wx wy wz ];
                    padding = 2*max( [ ws(~isinf(ws)) 1 ] );
                    if ~isinf(wx)
                        survey.int_xv = linspace( src.xrng(1)-padding, src.xrng(2)+padding );
                    else
                        survey.int_xv = 0;
                    end
                    if ~isinf(wy)
                        survey.int_yv = linspace( src.yrng(1)-padding, src.yrng(2)+padding );
                    else
                        survey.int_yv = 0;
                    end
                    if ~isinf(wz)
                        survey.int_zv = linspace( src.zrng(1)-padding, src.zrng(2)+padding );
                    else
                        survey.int_zv = 0;
                    end
                    [ survey.int_X, survey.int_Y, survey.int_Z ] = meshgrid( survey.xv, survey.yv, survey.zv );
                    survey.int_p = [ survey.X(:)'; survey.Y(:)'; survey.Z(:)' ];
                else
                    % Here, the user has directly supplied an array of points.
                    user_p = varargin{1};
                    dimsizes = size(user_p);
                    if dimsizes(1) ~= 3
                        error( 'Input list of evaluation points must be a 3xN array' );
                    else
                        survey.p = user_p;
                    end
                end
            end
        end

        % Setter functions.
        function thisSurvey = set.xv( thisSurvey, new_xv )
            % Updating xv means updating all downstream properties too.
            thisSurvey.int_xv = new_xv;
            [ thisSurvey.int_X, thisSurvey.int_Y, thisSurvey.int_Z ] = meshgrid( thisSurvey.xv, thisSurvey.yv, thisSurvey.zv );
            thisSurvey.int_p = [ thisSurvey.X(:)'; thisSurvey.Y(:)'; thisSurvey.Z(:)' ];
        end
        function thisSurvey = set.yv( thisSurvey, new_yv )
            % Updating yv means updating all downstream properties too.
            thisSurvey.int_yv = new_yv;
            [ thisSurvey.int_X, thisSurvey.int_Y, thisSurvey.int_Z ] = meshgrid( thisSurvey.xv, thisSurvey.yv, thisSurvey.zv );
            thisSurvey.int_p = [ thisSurvey.X(:)'; thisSurvey.Y(:)'; thisSurvey.Z(:)' ];
        end
        function thisSurvey = set.zv( thisSurvey, new_zv )
            % Updating zv means updating all downstream properties too.
            thisSurvey.int_zv = new_zv;
            [ thisSurvey.int_X, thisSurvey.int_Y, thisSurvey.int_Z ] = meshgrid( thisSurvey.xv, thisSurvey.yv, thisSurvey.zv );
            thisSurvey.int_p = [ thisSurvey.X(:)'; thisSurvey.Y(:)'; thisSurvey.Z(:)' ];
        end
        function thisSurvey = set.p( thisSurvey, new_p )
            thisSurvey.int_p = new_p;
            % If the user supplies p directly, then the xv, yv, zv, and
            % meshgrids are not relevant and should be blanked out to avoid
            % confusion.
            thisSurvey.int_xv = [];
            thisSurvey.int_yv = [];
            thisSurvey.int_zv = [];
            thisSurvey.int_X = [];
            thisSurvey.int_Y = [];
            thisSurvey.int_Z = [];
        end

        % Getter functions.
        function xv = get.xv( thisSurvey )
            xv = thisSurvey.int_xv;
        end
        function yv = get.yv( thisSurvey )
            yv = thisSurvey.int_yv;
        end
        function zv = get.zv( thisSurvey )
            zv = thisSurvey.int_zv;
        end
        function X = get.X( thisSurvey )
            X = thisSurvey.int_X;
        end
        function Y = get.Y( thisSurvey )
            Y = thisSurvey.int_Y;
        end
        function Z = get.Z( thisSurvey )
            Z = thisSurvey.int_Z;
        end
        function p = get.p( thisSurvey )
            p = thisSurvey.int_p;
        end
        function Np = get.Np( thisSurvey )
            Np = size( thisSurvey.p, 2 );
        end
        function dimsizes = get.dimsizes( thisSurvey )
            if ~isempty( thisSurvey.xv )
                dimsizes = [ length(thisSurvey.xv) length(thisSurvey.yv) length(thisSurvey.zv) ];
            end
        end
        function singdims = get.singdims( thisSurvey )
            singdims = thisSurvey.dimsizes == 1;
        end
        function nonsingdims = get.nonsingdims( thisSurvey )
            nonsingdims = thisSurvey.dimsizes > 1;
        end

        % Visualization.
        function overlaySurvey( thisSurvey, fh, varargin )
            % Overlay survey plane or transect for 2D or 1D surveys
            % (respectively). Do nothing for 3D surveys.
            assert( all(ishandle(fh)) );
            for k = 1 : length(fh)
                args = BaseTools.argarray2struct( [ { 'clr', [ .7 0 .7 ] } varargin ] );
                ah = fh(k).CurrentAxes;
                hold(ah,'on');
                if sum(thisSurvey.nonsingdims) == 1
                    plot3( ah, thisSurvey.p(1,:), thisSurvey.p(2,:), thisSurvey.p(3,:), '--', 'LineWidth', 2.0, 'Color', args.clr );
                elseif sum(thisSurvey.nonsingdims) == 2
                    sh = surf( ah, thisSurvey.X, thisSurvey.Y, thisSurvey.Z );
                    sh.EdgeColor = 'none';
                    sh.FaceColor = args.clr;
                    sh.FaceAlpha = 0.2;
                elseif sum(thisSurvey.nonsingdims) == 3
                    warning( 'Not implemented for 3D surveys.' );
                end
            end
        end

    end
end