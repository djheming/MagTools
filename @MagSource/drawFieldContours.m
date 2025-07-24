function ah = drawFieldContours( ah, survey, V, varargin )

% ah = drawFieldContours( ah, survey, V, ... ) make a figure showing
% contours of the supplied field.  
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
args = BaseTools.argarray2struct( [ { 'view', [], 'figsize', 520, 'overlay', true, 'cbar', true, 'suppress_2D_zero_contours', false, 'plane_shading_alpha', 0.7 } varargin ] );
if ~args.overlay
    cla( ah );
end

% Capture axis values and labels in cell arrays. This allows us to later
% index into these arrays flexibly because only one or two (not three) of
% them will be used for display.
axisvals = {
    squeeze(survey.xv);
    squeeze(survey.yv);
    squeeze(survey.zv); };
axislbls = {
    'x';
    'y'
    'z'
    };

% Need to reshape V.
V = reshape(V,size(survey.X));
V(isnan(V)) = 0; % Also, turn the NaNs into zeros.

% The rest depends on whether we're looking at a transect, a 2D plane, or a
% 3D volume (determined by the number of singeton vs nonsingleton dimensions).
if sum(survey.nonsingdims)==3

    % Here, we have a 3D volume. Use isosurfaces.
    def_view = [ 25 30 ];
    grid( ah, 'on' );
    axis( ah, 'equal' );
    xlabel( ah, 'x' );
    ylabel( ah, 'y' );
    zlabel( ah, 'z' );

    % Add isosurfaces for each contour value.
    n = 21;
    clrs = MagSource.magcolors(n);
    if isfield( args, 'ctrs' ) && ~isempty( args.ctrs )
        ctrs = args.ctrs;
    elseif isfield( args, 'clim' ) && ~isempty( args.clim )
        ctrs = asin(linspace(-1,1,n));
        ctrs = ctrs( ctrs>=min(args.clim) & ctrs<=max(args.clim) );
    elseif isfield( args, 'Vmax' ) && ~isempty( args.Vmax )
        ctrs = asin(linspace(-1,1,n)) * args.Vmax;
    else
        ctrs = asin(linspace(-1,1,n)) * max(abs(V(:)))/(2*pi);
    end
    for k = 1 : length(ctrs)
        s = isosurface( survey.X, survey.Y, survey.Z, V, ctrs(k) );
        ph = patch(s,'Parent',ah);
        ph.EdgeColor = 'none';
        if ctrs(k) == 0
            ph.FaceColor = [ .7 .7 .7 ];
        else
            ph.FaceColor = clrs(k,:);
        end
        ph.FaceAlpha = 0.25;
    end

    % Restrict axis limits to display range.
    ah.XLim = [ min(survey.xv) max(survey.xv) ];
    ah.YLim = [ min(survey.yv) max(survey.yv) ];
    ah.ZLim = [ min(survey.zv) max(survey.zv) ];

elseif sum(survey.nonsingdims)==2

    % Here we have a 2D plane. This is not compatible with
    % overlay mode, so we'll need to clear the axes.
    % TO DO: add the option to make this a true 2D plane filled with
    % imagesc rather than drawn using surfaces. 
    cla( ah );
    grid( ah, 'on' );
    hold( ah, 'on' );
    axis( ah, 'equal' );
    view( ah, 3 );
    xlabel( ah, 'x' );
    ylabel( ah, 'y' );
    zlabel( ah, 'z' );

    % We also need our arrays to be 2D.
    sqX = squeeze(survey.X);
    sqY = squeeze(survey.Y);
    sqZ = squeeze(survey.Z);
    sqV = squeeze(V);

    % Show colored surface.
    surface( sqX, sqY, sqZ, 'CData', sqV, 'FaceColor', 'texturemap', 'EdgeColor', 'none', 'FaceAlpha', 'texturemap', 'AlphaData', args.plane_shading_alpha, 'AlphaDataMapping', 'none', 'FaceLighting', 'none', 'EdgeLighting', 'none', 'AmbientStrength', 1 );
    colormap( ah, MagSource.magcolors );
    axis( ah, 'tight' );

    % Now we want to overlay contours. But which contour values?
    if isfield( args, 'clim' ) && ~isempty( args.clim )
        ctrvals = BaseTools.get_round_step_size( args.clim, 30, true );
    else
        ctrvals = BaseTools.get_round_step_size( V, 30, true );
    end
    if args.suppress_2D_zero_contours
        ctrvals = ctrvals(ctrvals~=0);
    end
    if args.cbar
        cbar = colorbar(ah);
        if isfield( args, 'cbar_label' )
            cbar.Label.String = args.cbar_label;
            % cbar.Label.Position = [ 0 1.1*max(ctrvals) 0 ];
            % cbar.Label.Rotation = 0;
            % cbar.Label.HorizontalAlignment = 'left';
        end
    end

    % Now overlay those contours. We need a custom function for this
    % because Matlab can only draw contours on the xy plane whereas, in
    % general, our contours could be on a plane perpendicular to any of the
    % three coordinate axes.
    if length(ctrvals)>1
        clim( ah, [ min(ctrvals) max(ctrvals) ] );
        plot_2D_contours_on_3D_plane( survey, V, ctrvals );
    else
        clim( ah, [ -1 1 ] );
    end

    % Set default view for 2D plane.
    switch axislbls{survey.singdims}
        case 'x'
            def_view = [ 90 0 ];
        case 'y'
            def_view = [ 0 0 ];
        case 'z'
            def_view = [ 0 90 ];
    end    

elseif sum(survey.nonsingdims)==1

    % Here, we will have a single line plot. Which dimension?
    switch axislbls{survey.nonsingdims}
        case 'x'
            ind = 1;
        case 'y'
            ind = 2;
        case 'z'
            ind = 3;
    end
    plot( ah, axisvals{ind}, V, 'LineWidth', 2.0 );
    xlabel( ah, axislbls{ind} );
    grid on;
    def_view = [ 0 90 ];
    
else
    error( 'Unexpected number of singleton dimensions.' );
end

% Adjust perspective, if necessary.
if isfield( args, 'view' ) && ~isempty( args.view )
    view( ah, args.view );
elseif exist( 'def_view', 'var' ) && ~isempty( def_view )
    view( ah, def_view );
else
    view( ah, [ 0 90 ] );
end

% If the user wishes a view where z points downward, we have to change
% the camera up property. This only works if the projection is not also
% orthographic.
% aligned = all(mod(ah.View,90)==0);
% if isfield( args, 'z_down' ) && args.z_down && ~aligned
%     camup( ah, [ 0 0 -1 ] );
% end
if isfield( args, 'z_down' ) && args.z_down
    ah.YDir = 'reverse';
    ah.ZDir = 'reverse';
end

% Illuminate the 3D scene.
camlight(ah);
lighting( ah, 'gouraud' );








function plot_2D_contours_on_3D_plane( survey, V, ctrvals )

% What is the orientation of our plane?
axislbls = {
    'x';
    'y'
    'z'
    };
switch axislbls{survey.singdims}
    case 'x'
        cV = squeeze(V)';
        cX = survey.yv;
        cY = survey.zv;
        X0 = survey.xv;
    case 'y'
        cV = squeeze(V)';
        cX = survey.xv;
        cY = survey.zv;
        Y0 = survey.yv;
    case 'z'
        cV = squeeze(V);
        cX = survey.xv;
        cY = survey.yv;
        Z0 = survey.zv;
end

% Make a contour matrix and extract the data.
Cmat = contourc( cX, cY, cV, ctrvals );
ctrs = BaseTools.parseContourMatrix( Cmat );

% Redraw each contour line, but now on the correct plane.
switch axislbls{survey.singdims}
    case 'x'
        for k = 1 : length(ctrs)
            plot3( X0*ones(size(ctrs(k).x)), ctrs(k).x, ctrs(k).y, 'k' );
        end
    case 'y'
        for k = 1 : length(ctrs)
            plot3( ctrs(k).x, Y0*ones(size(ctrs(k).x)), ctrs(k).y, 'k' );
        end
    case 'z'
        for k = 1 : length(ctrs)
            plot3( ctrs(k).x, ctrs(k).y, Z0*ones(size(ctrs(k).x)), 'k' );
        end
end




