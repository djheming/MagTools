function fh = showQfieldContours( thisSource, element, varargin )

% fh = showQfieldContours( thisSource, element, ... ) make a figure showing
% the specified element of the Q field as contours in a 3D volume.
%
%
%   Disclaimer: This code is provided as-is, has been tested only very
%   informally, and may not always behave as intended. I find it useful for
%   my own work, and I hope you will too, but I make no guarantees as to
%   the accuracy or robustness of this code. This code is also actively
%   under development and future versions may not be backward compatible.
%
%   Created: 2025-05-07
%   Doug Hemingway (douglas.hemingway@utexas.edu)
%   University of Texas at Austin
%


% Read input arguments.
args = BaseTools.argarray2struct( varargin, { 'view', [], 'figsize', 520, 'title', false, 'clim', [ -2*pi 2*pi ], 'show_source', true, 'show_2D_peaks_and_zeros', true } );
optargs = BaseTools.struct2argarray(args);

% Check leading positional arguments for a SurveyField. If it was not supplied, make it here.
for k = 1 : length(args.posArgs)
    if isa( args.posArgs{k}, 'SurveyField' )
        survey = args.posArgs{k};
    end
end
if ~exist( 'survey', 'var' ) || isempty( survey )
    survey = SurveyField( thisSource ); % Construct a survey for this source.
end

% If the user has supplied an array of surveys, handle them recursively.
num_surveys = length(survey);
if num_surveys > 1
    fh = cell(num_surveys,1);
    for k = 1 : num_surveys
        fh{k} = thisSource.showQfieldContours( element, survey(k), optargs{:} );
    end
    return;
end

% Compute Q fields (merging the contributions from each source).
Q = thisSource.computeQfield(survey.p,Qframe='A');

% Make shorthand for element letters.
eltxt = 'xyz';

% The user may also want to use the shorthand 'all' for all nine elements.
% We'll convert that internally to a more explicit list.
if ( ischar(element) || isstring(element) ) && strcmp( element, 'all' )
    element = {
        'xx' 'xy' 'xz' 
        'yx' 'yy' 'yz' 
        'zx' 'zy' 'zz' 
        };
end

% Capture axis values and labels in cell arrays. This allows us to later
% index into these arrays flexibly because some of them will be used for
% display. 
axisvals = {
    squeeze(survey.xv);
    squeeze(survey.yv);
    squeeze(survey.zv); };
axislbls = {
    'x';
    'y'
    'z'
    };

% If we have a 2D or 3D survey, we use drawFieldContours. If we have a 1D
% transect, we simply plot all the curves on the same axes.
if sum(survey.nonsingdims) > 1

    % Here, we have a 2D or 3D survey. And we may be asked to make plots
    % for any or all of the components.
    num_elements = numel(element);
    fh = gobjects(size(element));
    for k = 1 : num_elements

        % Establish the figure handle and grab the plotting axes.
        fh(k) = figure;
        ah = axes('Parent',fh(k));
        axis(ah,'equal'); % Since all dimensions of this plot are spatial.

        % Pull out the appropriate element and draw contours.
        i = eltxt_to_num( element{k}(1) );
        j = eltxt_to_num( element{k}(2) );
        MagSource.drawFieldContours( ah, survey, Q(i,j,:), optargs{:} );

        % Show the source body or bodies as well.
        if sum(survey.nonsingdims) > 1 && args.show_source
            thisSource.drawSource( ah, 'show_M', false );
            % If the source extends beyond the survey domain, it will have
            % expanded the figure. We'll need to scale back to the survey
            % domain.
            if survey.nonsingdims(1)
                ah.XLim = [ min(survey.xv) max(survey.xv) ];
            end
            if survey.nonsingdims(2)
                ah.YLim = [ min(survey.yv) max(survey.yv) ];
            end
            if survey.nonsingdims(3)
                ah.ZLim = [ min(survey.zv) max(survey.zv) ];
            end
        end

        % Optionally add a title.
        if args.title
            text( 0, max(survey.Y(:)), max(survey.Z(:)), [ 'Q_{' eltxt(i) eltxt(j) '}' ], 'FontSize', 16 );
        end

        % If this is a 2D survey for a prism, there are some optional overlays.
        if sum(survey.nonsingdims) == 2 && isa(thisSource,'MagBox') && sum(thisSource.finitedims)==2

            % Optionally add panel label.
            if isfield( args, 'axlbl' ) && ~isempty( args.axlbl )
                text( 0.03, 0.97, args.axlbl, 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 16 );
            end

            % Optionally show zero crossings and peaks.
            if args.show_2D_peaks_and_zeros
                switch element{k}
                    % TO DO: right now, I pass in survey.xv but this may
                    % not always be right. Think of a more general
                    % solution.
                    case { 'xx', 'yy', 'zz' }
                        pQiizc = thisSource.computeQiiZeros( survey.xv );
                        plot3( ah, pQiizc(1,:), pQiizc(2,:), pQiizc(3,:), 'k--', 'LineWidth', 1.5 );
                    otherwise
                        pQijzc = thisSource.computeQijZeros( survey.xv );
                        plot3( ah, pQijzc(1,:), pQijzc(2,:), pQijzc(3,:), 'k:', 'LineWidth', 1.5 );
                        [ p_short, p_long ] = thisSource.computeQijPeaks( survey.xv );
                        plot3( ah, p_short(1,:), p_short(2,:), p_short(3,:), 'Color', [ .4 .4 .7 ], 'LineStyle', '-.', 'LineWidth', 1.5 );
                        plot3( ah, p_long(1,:), p_long(2,:), p_long(3,:), 'Color', [ .7 .4 .4 ], 'LineStyle', '-.', 'LineWidth', 1.5 );
                end
            end
            
        end
        
        % Resize.
        fh(k).Position(3:4) = args.figsize;

    end

elseif sum(survey.nonsingdims) == 1
    
    % Here, we have a 1D survey. Draw a line for each element, all on the
    % same axes.
    num_elements = numel(element);
    clrs = lines(3);
    styls = { '-', '-.', '--' };
    hvals = axisvals{survey.nonsingdims};
    fh = figure;
    ah = axes('Parent',fh);
    hold(ah,'on');
    grid(ah,'on');
    legtxt = cell(num_elements,1);
    element = element'; % Force row-ordering.
    for k = 1 : num_elements
        i = eltxt_to_num( element{k}(1) );
        j = eltxt_to_num( element{k}(2) );
        legtxt{k} = [ 'Q_{' eltxt(i) eltxt(j) '}' ];
        plot( ah, hvals, squeeze(Q(i,j,:)), 'LineStyle', styls{j}, 'Color', clrs(i,:), 'LineWidth', 2.0 );
    end
    xlabel( axislbls{survey.nonsingdims} );
    ylabel( 'B (nT)' );
    lh = legend( legtxt );
    lh.AutoUpdate = 'off';
    
else
    warning( 'Cannot show Q-field contours without a regular survey grid.' );
end    


function i = eltxt_to_num( txt )

switch txt
    case 'x'
        i = 1;
    case 'y'
        i = 2;
    case 'z'
        i = 3;
    otherwise
        i = 0;
end




