function ah = drawBox( thisMagBox, varargin )

% ah = drawBox( thisMagBox, ... ) draw shaded rectangles showing
% the faces of the magnetized prism represented by thisMagBox.
%
%
%   Disclaimer: This code is provided as-is, has been tested only very
%   informally, and may not always behave as intended. I find it useful for
%   my own work, and I hope you will too, but I make no guarantees as to
%   the accuracy or robustness of this code. This code is also actively
%   under development and future versions may not be backward compatible.
%
%   Created: 2025-04-14
%   Doug Hemingway (douglas.hemingway@utexas.edu)
%   University of Texas at Austin
%


% Read input arguments.
def_args = { 'overlay', true, 'axislabels', 'xyz', 'box_clr', [ .4 .4 .4 ], 'axes', false, 'label_origin', false, 'show_p_vector', false, 'label_p_vector', false, 'label_p', false, 'faces', true, 'vertices', false, 'verts_from_origin', [], 'rlabels', true, 'show_M', false, 'M_length', [] };
args = BaseTools.argarray2struct( varargin, def_args );

% Check leading positional arguments for an axis handle.
for k = 1 : length(args.posArgs)
    if isa( args.posArgs{k}, 'matlab.graphics.axis.Axes' )
        ah = args.posArgs{k};
    end
end

% Other settings.
p_clr = [ 0 .4 0 ];
pS_clr = [ 0 .6 0 ];
FS_clr = [ .4 .4 .4 ];

% Verify that we have a plot axis and bring it into focus.
if ~exist( 'ah', 'var' ) || ~ishandle( ah )
    fh = figure;
    ah = axes('Parent',fh);
    xlabel(args.axislabels(1));
    ylabel(args.axislabels(2));
    zlabel(args.axislabels(3));
else
    % Remember the axis limits and go back to them at the end.
    axlims = axis(ah);
end
hold(ah,'on');
grid(ah,'on');
axis(ah,'tight');
axis(ah,'equal');
if ~args.overlay
    cla( ah );
    colorbar( ah, 'off' );
end

% If the user wishes to show the coordinate frames, do that here.
if args.axes

    % Show analysis coordinate frame.
    BaseTools.drawFrame( ah, zeros(3,1), eye(3), 'k', args.axislabels );
    if args.label_origin
        text( 0, -.2, 0, '(A)', 'FontSize', 14 );
    end

    % If there is a distinct source coordinate frame, show that too.
    if ( ~isempty(thisMagBox.v_AS) && ~all(thisMagBox.v_AS==zeros(3,1)) ) || ~all(all((thisMagBox.R_AS==eye(3))))

        % Note that we want to feed this function R_SA (not R_AS)
        % because we're not trying to transform from (S) to (A),
        % we're only trying to display the components of the (S)
        % frame.
        if args.axes
            BaseTools.drawFrame( ah, thisMagBox.v_AS, thisMagBox.R_SA, FS_clr, args.axislabels );
            if args.label_origin
                text( ah, thisMagBox.v_AS(1), thisMagBox.v_AS(2)-.2, thisMagBox.v_AS(3), '(S)', 'FontSize', 14, 'Color', FS_clr );
            end
        end

        % Draw vector connecting the two frames.
        BaseTools.drawArrow( ah, [ 0 thisMagBox.v_AS(1) ], [ 0 thisMagBox.v_AS(2) ], [ 0 thisMagBox.v_AS(3) ], 'Color', FS_clr, 'LineWidth', 2.0, 'LineStyle', '--', 'headlength', 0.1 );
        if args.label_origin
            text( ah, 0.6*thisMagBox.v_AS(1), 0.6*thisMagBox.v_AS(2), 0.6*thisMagBox.v_AS(3), 'v^{AS}', 'Color', FS_clr, 'FontWeight', 'bold', 'FontSize', 14, 'HorizontalAlignment', 'right' );
        end

    end

end

% If the source is infinite in one of its directions, we'll artificially
% flatten it first.
if any(isinf(abs(thisMagBox.vx)))
    tmprng.x = 0;
else
    tmprng.x = thisMagBox.vx;
end
if any(isinf(abs(thisMagBox.vy)))
    tmprng.y = 0;
else
    tmprng.y = thisMagBox.vy;
end
if any(isinf(abs(thisMagBox.vz)))
    tmprng.z = 0;
else
    tmprng.z = thisMagBox.vz;
end

% If the source is finite in volume, it's easy enough to draw its faces.
% But don't forget to transform from the prism coordinates to the analysis
% coordinate system.
if sum(thisMagBox.finitedims)==3

    % X-faces. 
    for j = 1 : 2
        [ Xs, Ys, Zs ] = meshgrid( thisMagBox.vx(j), thisMagBox.vy, thisMagBox.vz );
        [ X, Y, Z ] = transform_coordinates( thisMagBox, Xs, Ys, Zs );
        sh = surf( ah, squeeze(X), squeeze(Y), squeeze(Z) );
        sh.EdgeColor = args.box_clr;
        sh.LineWidth = 1.0;
        if args.faces
            sh.FaceColor = args.box_clr;
        else
            sh.FaceColor = 'none';
        end
        sh.FaceAlpha = 0.3;
    end

    % Y-faces.
    for j = 1 : 2
        [ Xs, Ys, Zs ] = meshgrid( thisMagBox.vx, thisMagBox.vy(j), thisMagBox.vz );
        [ X, Y, Z ] = transform_coordinates( thisMagBox, Xs, Ys, Zs );
        sh = surf( ah, squeeze(X), squeeze(Y), squeeze(Z) );
        sh.EdgeColor = args.box_clr;
        sh.LineWidth = 1.0;
        if args.faces
            sh.FaceColor = args.box_clr;
        else
            sh.FaceColor = 'none';
        end
        sh.FaceAlpha = 0.3;
    end

    % Z-faces.
    for j = 1 : 2
        [ Xs, Ys, Zs ] = meshgrid( thisMagBox.vx, thisMagBox.vy, thisMagBox.vz(j) );
        [ X, Y, Z ] = transform_coordinates( thisMagBox, Xs, Ys, Zs );
        sh = surf( ah, squeeze(X), squeeze(Y), squeeze(Z) );
        sh.EdgeColor = args.box_clr;
        sh.LineWidth = 1.0;
        if args.faces
            sh.FaceColor = args.box_clr;
        else
            sh.FaceColor = 'none';
        end
        sh.FaceAlpha = 0.3;
    end

elseif sum(thisMagBox.finitedims)==2

    % If the sources is infinitely long in one dimension, we will not be able
    % to draw its six faces in the usual way so we'll need to truncate.
    [ Xs, Ys, Zs ] = meshgrid( tmprng.x, tmprng.y, tmprng.z );
    [ X, Y, Z ] = transform_coordinates( thisMagBox, Xs, Ys, Zs );
    sh = surf( ah, squeeze(X), squeeze(Y), squeeze(Z) );
    sh.EdgeColor = args.box_clr;
    sh.LineWidth = 1.0;
    if args.faces
        sh.FaceColor = args.box_clr;
    else
        sh.FaceColor = 'none';
    end
    sh.FaceAlpha = 0.3;

end
if exist( 'axlims', 'var' ) && ~isempty( axlims )
    axis( ah, axlims );
end

% Does the user have a preferred viewing orientation?
if isfield( args, 'view' ) && ~isempty( args.view )
    view( ah, args.view );
    camlight; % Give it some 3D lighting.
end

% Does the user wish to display the magnetization direction?
if args.show_M
    lenM = norm(thisMagBox.M);
    if ~isempty( args.M_length )
        lenM = lenM/args.M_length;
    end
    BaseTools.drawArrow( ah, thisMagBox.x0+[ 0 thisMagBox.M(1)/lenM ], thisMagBox.y0+[ 0 thisMagBox.M(2)/lenM ], thisMagBox.z0+[ 0 thisMagBox.M(3)/lenM ], 'Color', 'k', 'LineWidth', 2.0 );
end

% If the user wishes to show vectors connecting the origin to some of the
% vertices, do that here.
if isfield( args, 'verts_from_origin' ) && ~isempty( args.verts_from_origin ) && iscell( args.verts_from_origin )
    v_clr = [ .4 0 .4 ];
    for k = 1 : length(args.verts_from_origin)
        a = str2double(args.verts_from_origin{k}(1));
        b = str2double(args.verts_from_origin{k}(2));
        c = str2double(args.verts_from_origin{k}(3));
        BaseTools.drawArrow( ah, [ 0 thisMagBox.xrng(a) ], [ 0 thisMagBox.yrng(b) ], [ 0 thisMagBox.zrng(c) ], 'Color', v_clr, 'LineWidth', 2.0, 'headlength', 0.1 );
    end
end

% If the user has also provided an evaluation point, display that too.
if isfield( args, 'p' )

    % Check size of p.
    [ rows, cols ] = size(args.p);
    if rows ~= 3
        error( 'p must have 3 rows' );
    end
    if cols > 1
        error( 'Cannot show more than 1 p vector at a time' );
    end

    % Draw p vector from analysis coordinate system origin.
    if args.show_p_vector
        BaseTools.drawArrow( ah, [ 0 args.p(1) ], [ 0 args.p(2) ], [ 0 args.p(3) ], 'Color', p_clr, 'LineWidth', 2.0, 'headlength', 0.1 );
        if isfield( args, 'label_p_vector' ) && args.label_p_vector
            text( ah, 0.8*args.p(1), 0.8*args.p(2), 0.8*args.p(3), 'p^{(A)}', 'Color', p_clr, 'FontWeight', 'bold', 'FontSize', 14, 'HorizontalAlignment', 'right' );
        elseif args.label_p
            text( 1.05*args.p(1), 1.05*args.p(2), 1.05*args.p(3), 'p', 'FontSize', 14, 'FontWeight', 'bold', 'Color', p_clr, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 14 );
        end
    end

    % Might also need to draw the vector from the source coordinate frame.
    if ~isempty(thisMagBox.v_AS) && ~all(thisMagBox.v_AS==zeros(3,1))
        BaseTools.drawArrow( ah, [ thisMagBox.v_AS(1) args.p(1) ], [ thisMagBox.v_AS(2) args.p(2) ], [ thisMagBox.v_AS(3) args.p(3) ], 'Color', pS_clr, 'LineWidth', 2.0, 'headlength', 0.1 );
        if isfield( args, 'label_p_vector' ) && args.label_p_vector
            pS = args.p - thisMagBox.v_AS;
            text( ah, 0.8*pS(1)+thisMagBox.v_AS(1), 0.8*pS(2)+thisMagBox.v_AS(2), 0.8*pS(3)+thisMagBox.v_AS(3), 'p^{(S)}', 'Color', pS_clr, 'FontWeight', 'bold', 'FontSize', 14, 'HorizontalAlignment', 'right' );
        end
    end

    % Draw lines connecting each vertex to p.
    if isfield( args, 'vertices' ) && args.vertices

        % Obtain the vertex coordinates in the Analysis frame. If any of
        % the vertices are at infinity, we won't be able to draw that so
        % let's flatten those to zeros. In the case that the prism has
        % infinite extent in one of its dimensions, we can flatten that
        % dimension so there are effectively only 4 vertices rather than 8.
        % We accomplish this by making the ranges for each of a,b,c depend
        % on whether or not the corresponding dimension is finite.
        v = reshape( thisMagBox.vA, 3, 2, 2, 2 );
        v(isinf(v)) = 0;
        for a = unique( [ 1 1 + thisMagBox.finitedims(1) ] )
            rclr = [ 0 0.25*(1+a) 1 ];
            for b = unique( [ 1 1 + thisMagBox.finitedims(2) ] )
                for c = unique( [ 1 1 + thisMagBox.finitedims(3) ] )
                    line( ah, [ v(1,a,b,c) args.p(1) ], [ v(2,a,b,c) args.p(2) ], [ v(3,a,b,c) args.p(3) ], 'LineWidth', 2.0, 'Color', rclr );
                    if isfield( args, 'rlabels' ) && args.rlabels
                        abcstr = [ num2str(a) num2str(b) num2str(c) ]; 
                        rlbl = [ 'r_{' abcstr(thisMagBox.finitedims) '}' ]; % Select only the characters cooresponding to finite dimensions.
                        text( ah, v(1,a,b,c), v(2,a,b,c), v(3,a,b,c), rlbl, 'Color', rclr, 'FontWeight', 'bold', 'FontSize', 14 );
                    end
                end
            end
        end

        % If p is not already shown, show it here in the same color as the
        % r vectors.
        if ~args.show_p_vector
            text( args.p(1), args.p(2), args.p(3), 'p', 'FontSize', 14, 'FontWeight', 'bold', 'Color', [ 0 0.5 1 ], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14 );
        end

        % Might also want to shade the resulting pyramids?
        % Note: For now, this is implemented only for the ±x faces and only
        % works when the prism has finite volume.
        if isfield( args, 'shading' ) && args.shading
            if sum(thisMagBox.finitedims)==3
                for a = 1 : 2
                    fclr = [ 0 0.25*(1+a) 1 ];
                    for b = 1 : 2
                        for c = 1 : 2
                            % If odd, flip c; if even, flip b.
                            if mod(sum( [ a b c ] ),2)
                                b2 = b;
                                c2 = mod(c,2)+1;
                            else
                                b2 = mod(b,2)+1;
                                c2 = c;
                            end
                            X = [ v(1,a,b,c) v(1,a,b2,c2) args.p(1) ]';
                            Y = [ v(2,a,b,c) v(2,a,b2,c2) args.p(2) ]';
                            Z = [ v(3,a,b,c) v(3,a,b2,c2) args.p(3) ]';
                            ph = patch( ah, X, Y, Z, fclr );
                            ph.EdgeColor = 'none';
                            ph.FaceAlpha = 0.3;
                        end
                    end
                end
            else
                error( 'Not implemented for prisms with infinite extent.' );
            end
        end

        % Might also want to show the alphas and f ratios?
        % Note: For now, this is implemented only for i=x, j=y.
        if isfield( args, 'alphas_and_fs' ) && args.alphas_and_fs
            
            af_txt_clr = [ 0 0.5 1 ];

            % Compute r vectors.
            rx1 = args.p(1) - min(thisMagBox.vx);
            rx2 = args.p(1) - max(thisMagBox.vx);
            ry1 = args.p(2) - min(thisMagBox.vy);
            ry2 = args.p(2) - max(thisMagBox.vy);
            rz1 = args.p(3) - min(thisMagBox.vz);
            rz2 = args.p(3) - max(thisMagBox.vz);

            % Compute Qii and Qij along with alphas and f ratios.
            [ Qxx, alphas, ~ ] = MagBox.computeQii( [ rx1; rx2 ], [ ry1; ry2 ], [ rz1; rz2 ] );
            [ Qxy, f_ab, ~, ~, Rsq ] = MagBox.computeQij( [ rx1; rx2 ], [ ry1; ry2 ], [ rz1; rz2 ] );

            % Add text to show the relevant alpha and length ratio values.
            % The way we do this depends on whether or not the prism is
            % finite. If the prism has infinite extent in the k-direction,
            % we have to modify things slightly because the alphas have a
            % slightly different definition (they're half of what they
            % are in the finite volume case) and instead of f ratios, we go
            % directly to the r_{ab} lengths.
            if sum(thisMagBox.finitedims)==3
                Qii_txt = [ '\alpha_{1}=' num2str(alphas(1),3) ', \alpha_{2}=' num2str(alphas(2),3) ', Q_{ii}=' num2str(Qxx,3) ];
                Qij_txt = [ 'f_{1,1}=' num2str(f_ab(1,1),3) ', f_{1,2}=' num2str(f_ab(1,2),3) ', f_{2,1}=' num2str(f_ab(2,1),3) ', f_{2,2}=' num2str(f_ab(2,2),3) ', Q_{ij}=' num2str(Qxy,3) ];
            elseif ~thisMagBox.finitedims(3)
                alphas = alphas/2;
                Qii_txt = [ '\alpha_{1}=' num2str(alphas(1),3) ', \alpha_{2}=' num2str(alphas(2),3) ', Q_{ii}=' num2str(Qxx,3) ];
                Qij_txt = [ '|r_{1,1}|=' num2str(sqrt(Rsq(1,1)),3) ', |r_{1,2}|=' num2str(sqrt(Rsq(1,2)),3) ', |r_{2,1}|=' num2str(sqrt(Rsq(2,1)),3) ', |r_{2,2}|=' num2str(sqrt(Rsq(2,2)),3) ', Q_{ij}=' num2str(Qxy,3) ];
                % In the 2D case, we can also show shaded angles for the alphas.
                theta11 = atan2( ry1, rx1 );
                theta12 = atan2( ry2, rx1 );
                theta21 = atan2( ry1, rx2 );
                theta22 = atan2( ry2, rx2 );
                overlay_shaded_angle( theta11, theta12, args.p(1), args.p(2), 1, [ 0 0.5 1 ], '\alpha_1' );
                overlay_shaded_angle( theta22, theta21, args.p(1), args.p(2), .8, [ 0 0.5 1 ], '\alpha_2' );
                delete(findall(gcf,'Type','light'));
            else
                error( 'Not implemented for infinite extent in x or y dimensions.' );
            end
            h = annotation( 'textbox', [0, 0, 0.1, 0.1], ...
                'String', { Qii_txt, Qij_txt }, ...
                'Color', af_txt_clr, 'HorizontalAlignment', 'center', ...
                'FitBoxToText', 'on', 'BackgroundColor', 'none', 'EdgeColor', 'none' );
            drawnow;
            h.Position(1) = 0.5 - h.Position(3)/2;
            h.Position(2) = 0.8;
            
        end

        % For the 2D case, the user may also want to show the zero
        % crossings.
        if isfield( args, 'zerocrossings' ) && args.zerocrossings && sum(thisMagBox.finitedims)==2
            if ~thisMagBox.finitedims(3)
                b = max(thisMagBox.w(thisMagBox.finitedims));
                pQiizc = thisMagBox.computeQiiZeros( linspace(-b,b,41) );
                plot( ah, pQiizc(1,:), pQiizc(2,:), 'k--', 'LineWidth', 1.5 );
                pQijzc = thisMagBox.computeQijZeros( linspace(-b,b,41) );
                plot( ah, pQijzc(1,:), pQijzc(2,:), 'k:', 'LineWidth', 1.5 );
            else
                error( 'Not implemented for infinite extent in x or y dimensions.' );
            end
        end

    end

end




function [ Xa, Ya, Za ] = transform_coordinates( thisMagBox, Xs, Ys, Zs )

dimsizes = size(Xs);
pS = [ Xs(:)'; Ys(:)'; Zs(:)' ];
pA = thisMagBox.v_AS + thisMagBox.R_AS * pS;
Xa = reshape(pA(1,:),dimsizes);
Ya = reshape(pA(2,:),dimsizes);
Za = reshape(pA(3,:),dimsizes);




function overlay_shaded_angle( theta1, theta2, x0, y0, r, clr, lbl )

% Adjust angles for plotting: theta1 is the starting angle, theta2 is the
% ending angle. In general, we want to plot the opposite angle, so we need
% to add pi to both angles to get the correct azimuths for plotting.
theta1 = theta1 + pi;
theta2 = theta2 + pi;

% x0, y0 define the center of the circle
% r is the radius of the arc
theta = linspace( theta1, theta2 );
ori = sign(theta2-theta1);
x = x0 + r * cos(theta);
y = y0 + r * sin(theta);
BaseTools.drawArrow( gca, x, y, 'Color', clr, 'headlength', r/8 );

% Also draw as a shaded patch.
ph = fill( [ x x0 ], [ y y0 ], 'k' );
ph.EdgeColor = 'none';
ph.FaceColor = clr;
ph.FaceAlpha = 0.4 + 0.2*ori;

% Label (just beyond the range).
theta3 = theta2 - ori*0.15;
xe = x0 + r*1.3 * cos(theta3);
ye = y0 + r*1.3 * sin(theta3);
text( xe, ye, lbl, 'Color', clr, 'FontSize', 14, 'HorizontalAlignment', 'center' );
