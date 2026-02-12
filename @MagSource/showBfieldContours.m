function fh = showBfieldContours( thisSource, component, varargin )

% fh = showBfieldContours( thisSource, component, ... ) make a
% figure showing the specified element of the B field as contours in a 3D volume.
%
%
%   Disclaimer: This code is provided as-is, has been tested only very
%   informally, and may not always behave as intended. I find it useful for
%   my own work, and I hope you will too, but I make no guarantees as to
%   the accuracy or robustness of this code. This code is also actively
%   under development and future versions may not be backward compatible.
%
%   Created: 2025-05-09
%   Doug Hemingway (douglas.hemingway@utexas.edu)
%   University of Texas at Austin
%


% Read input arguments.
[ args, optargs ] = BaseTools.argarray2struct( varargin, { 'Q', [], 'view', [], 'z_down', false, 'figsize', 520, 'title', false, 'overlay_vectors', true, 'show_source', true } );

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
        fh{k} = thisSource.showBfieldContours( component, survey(k), optargs{:} );
    end
    return;
end

% Compute B field at all evaluation points.
B = thisSource.computeBfield(survey.p,'Q',args.Q);

% The intensity of magnetization is going to affect the contours we're
% going to want to plot. This is just a rule of thumb and can be overridden
% by having the calling function set Vmax.
if ~isfield( args, 'Vmax' ) || isempty( args.Vmax )
    args.Vmax = 100 * norm(thisSource.M);
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
    num_comps = numel(component);
    fh = gobjects(1,num_comps);
    for k = 1 : num_comps

        % Establish the figure handle and grab the plotting axes.
        fh(k) = figure;
        ah = axes('Parent',fh(k));
        axis(ah,'equal'); % Since all dimensions of this plot are spatial.

        % Pull out the appropriate B component and draw contours.
        Bcomp = thisSource.getBcomponent( B, component(k) );
        MagSource.drawFieldContours( ah, survey, Bcomp*1e9, 'cbar_label', [ 'B_{' component(k) '} (nT)' ], 'view', args.view, optargs{:}, 'Vmax', args.Vmax );

        % Show the source body or bodies as well.
        if sum(survey.nonsingdims) > 1 && args.show_source
            thisSource.drawSource( ah, optargs{:} );
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
            xtxt = max(survey.X(:));
            ytxt = max(survey.Y(:));
            if args.z_down
                ztxt = min(survey.Z(:));
            else
                ztxt = max(survey.Z(:));
            end
            text( 0, ytxt, ztxt, [ 'B_{' component(k) '}' ], 'FontSize', 16 );
        end

        % If this is a 2D survey, there are some optional overlays.
        if sum(survey.nonsingdims) == 2

            % Optionally overlay vector field.
            if isfield( args, 'overlay_vectors' ) && ~isempty( args.overlay_vectors ) && args.overlay_vectors
                % First, scale down the grid resolution to something
                % reasonable.
                xv = BaseTools.get_round_step_size( survey.xv, 16 );
                yv = BaseTools.get_round_step_size( survey.yv, 16 );
                zv = BaseTools.get_round_step_size( survey.zv, 16 );
                vector_survey = SurveyField( xv, yv, zv );
                Bv = thisSource.computeBfield(vector_survey.p);
                MagSource.drawBfieldVectors( ah, vector_survey, Bv, optargs{:} );
            end

            % Optionally add panel label.
            if isfield( args, 'axlbl' ) && ~isempty( args.axlbl )
                text( 0.03, 0.97, args.axlbl, 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 16 );
            end

        end
        
        % Resize.
        fh(k).Position(3:4) = args.figsize;

    end
    
elseif sum(survey.nonsingdims) == 1
    
    % How many components? Set line colors accoringly.
    num_comps = length(component);
    clrs = lines(num_comps);

    % Here, we have a 1D survey. Draw a line for each component, all on the
    % same axes.
    hvals = axisvals{survey.nonsingdims};
    fh = figure;
    ah = axes('Parent',fh);
    hold(ah,'on');
    grid(ah,'on');
    legtxt = cell(num_comps,1);
    for k = 1 : num_comps
        i = eltxt_to_num( component(k) );
        legtxt{k} = [ 'B_' component(k) ];
        plot( ah, hvals, squeeze(B(i,:))*1e9, 'Color', clrs(k,:), 'LineWidth', 2.0 );
    end
    xlabel( axislbls{survey.nonsingdims} );
    ylabel( 'B (nT)' );
    lh = legend( legtxt );
    lh.AutoUpdate = 'off';

    % Optionally annotate the 1D transect lines.
    if isfield( args, 'transect_labels' ) && ~isempty( args.transect_labels ) && length(args.transect_labels)==num_comps
        for k = 1 : num_comps
            % Add label at an approximately middle value of the curve.
            i = eltxt_to_num( component(k) );
            ev = abs(B(i,:)-mean(B(i,:)));
            [ ~, ind ] = min(ev);
            text( ah, hvals(ind), B(i,ind)*1e9, [ ' ' args.transect_labels{k} ], 'FontSize', 18, 'Color', clrs(k,:) );
        end
    end

    % Optionally overlay the field from a reference source.
    if isfield( args, 'ref_src' ) && ~isempty( args.ref_src ) && isa( args.ref_src, 'MagSource' )
        Bref = args.ref_src.computeBfield(survey.p);
        for k = 1 : num_comps
            i = eltxt_to_num( component(k) );
            plot( ah, hvals, squeeze(Bref(i,:))*1e9, 'Color', clrs(k,:), 'LineStyle', '--', 'LineWidth', 1.0 );
        end
    end

else
    warning( 'Cannot show B-field contours without a regular survey grid.' );
end


function i = eltxt_to_num( txt )

switch txt
    case 'x'
        i = 1;
    case 'y'
        i = 2;
    case 'z'
        i = 3;
    case 'xy'
        i = 1:2;
    case 'yz'
        i = 2:3;
    case 'xz'
        i = [ 1 3 ];
    case 't'
        i = 1:3;
    otherwise
        i = 0;
end

