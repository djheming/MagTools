function make_figures_for_paper( figs_folder )

% Make figures for the paper accompanying this software.
%
%
%   Disclaimer: This code is provided as-is, has been tested only very
%   informally, and may not always behave as intended. I find it useful for
%   my own work, and I hope you will too, but I make no guarantees as to
%   the accuracy or robustness of this code. This code is also actively
%   under development and future versions may not be backward compatible.
%
%   Doug Hemingway (douglas.hemingway@utexas.edu)
%   University of Texas at Austin
%   2024-06-15
%

% Make sure we have a destination for the output figures.
if ~exist( 'figs_folder', 'var' ) || isempty( figs_folder )
    figs_folder = '../Figures';
end
if ~exist( figs_folder, 'dir' )
    error( 'Please specify a valid destination folder for output figures.' );
end

% Select which plots to generate by setting the appropriate variables to true.
show_2D_box_geometry = false;
show_3D_box_geometry = false;
show_3D_coordinate_change = false;
show_2D_Q_and_B = true;
show_3D_Q_and_B = false;
show_Blakely = false;
show_Bongiolo = false;
show_wire_plate_limits = false;
show_wire_plate_movie = false;
show_shape_movie = false;
show_direction_movie = false;
show_semi_infinite_movies = false;

% Generate selected plots.
close all;
if show_2D_box_geometry
    
    % Create a 2D prism (infinitely long in the other dimension.
    wy = 2;
    wx = 4;
    myBox = MagBox( [ -wx/2 wx/2 ], [ -wy/2 wy/2 ], [ -Inf Inf ], [ 1 0 0 ]' );

    % Consider various evaluation points.
    p = [
        3 2 0
        -0.5 2.5 0
        ]';
    N = size(p,2);

    % Show the box in 2D, illustrating the evaluation point.
    for k = 1 : N
        ah = myBox.drawBox( 'p', p(:,k), 'label_p', true, 'vertices', true, 'alphas_and_fs', true, 'zerocrossings', true, 'view', [ 0 90 ] );
        exportgraphics( ancestor(ah,'figure'), sprintf('../Figures/2D_box_geometry_p(%s,%s).png', num2str(p(1,k)), num2str(p(2,k)) ) );
    end

end
if show_3D_box_geometry

    %
    % First, a simple illustration showing the r vectors.
    %

    % Establish a magnetized box and define an example evaluation point.
    myBox = MagBox( [ 1.5 3.5 ], [ 0.5 1.5 ], [ 0.3 1.3 ], [ 1 0 0 ]' );
    p = [ 2.45 -1 0.9 ]';

    % Show the box in 3D, illustrating the evaluation point.
    ah = myBox.drawBox( 'p', p, 'axes', true, 'vertices', true, 'show_p_vector', true, 'verts_from_origin', { '122', '211' }, 'view', [ 39 47 ] );

    % Save orthographic and oblique views.
    view( ah, [ 0 90 ] ); % xy
    exportgraphics( ah.Parent, '../Figures/3D_box_geometry_xy.png', 'ContentType', 'image', 'Resolution', 300 );
    view( ah, [ 0 0 ] ); % xz
    exportgraphics( ah.Parent, '../Figures/3D_box_geometry_xz.png', 'ContentType', 'image', 'Resolution', 300 );
    view( ah, [ 90 0 ] ); % yz
    exportgraphics( ah.Parent, '../Figures/3D_box_geometry_yz.png', 'ContentType', 'image', 'Resolution', 300 );
    view( ah, [ 39 47 ] );
    exportgraphics( ah.Parent, '../Figures/3D_box_geometry_oblique.png', 'ContentType', 'image', 'Resolution', 300 );


    %
    % Now show solid angles and line lengths to show geometric
    % representation of how Qii and Qij are determined.
    %

    % Make a box and define a couple of different evaluation points.
    myBox = MagBox( [ 1.5 3.5 ], [ 0.5 1.5 ], [ -1 0 ], [ 1 0 0 ]' );
    p = [
        2.45 -1 0.5
        4.2 -1 -0.4
        ]';
    num_points = size(p,2);

    % For each evaluation point, show the box with shaded solid angles on
    % the +i and -i faces.
    for k = 1 : num_points
        myBox.drawBox( 'p', p(:,k), 'vertices', true, 'show_p_vector', true, 'shading', true, 'axes', true, 'axislabels', 'ijk', 'rlabels', false, 'alphas_and_fs', true, 'view', [ 20 60 ] );
        axis( [ 0 4.5 -1 2 -1.5 2 ] );
        exportgraphics( gcf, [ '../Figures/box_geometry_with_lines_and_angles_' num2str(k) '.png' ] );
    end

end
if show_3D_coordinate_change

    % Define a prism that is a tad askew.
    roll = 10;
    pitch = -20;
    yaw = 30;
    R_AS = BaseTools.rpy2rot( roll, pitch, yaw );
    v_AS = [ 2 -.5 0 ]';
    M = [ 0 1 0 ]'; % This is in the source coordinate frame.
    myBox = MagBox( [ 0 2 ], [ 0 1 ], [ 0 0.6 ], M, v_AS, R_AS );

    % Show box with evaluation point.
    p = [ 2.5 1.5 0 ]';
    ah = myBox.drawBox( 'axes', true, 'p', p, 'label_origin', true, 'show_p_vector', true, 'label_p_vector', true, 'view', [ 0 90 ] );
    axis( ah, 'off' );
    camlight;
    exportgraphics( gcf, '../Figures/3D_coordinate_change.png', 'ContentType', 'image', 'Resolution', 300 );

end
if show_2D_Q_and_B

    % For the 2D problem, make a box that has infinite extent in the
    % z-direction.
    myBox = MagBox( [ -3 3 ], [ -1 1 ], [ -Inf Inf ], [ 3 1 0 ]' );
    fine_2D_survey = SurveyField( linspace(-12,12), linspace(-12,12), 0 );

    % Show Q field.
    fh = BaseTools.tile_figures( myBox.showQfieldContours( { 'xx', 'xy'; 'yx' 'yy' }, fine_2D_survey, 'title', true ) );
    fh.Position(4) = 920;
    exportgraphics( fh, '../Figures/2D_Q.png', 'ContentType', 'image', 'Resolution', 300 );

    % Now show B field vectors over a grid, but highlighting the vectors
    % along the Qii and Qij zero crossings. Let's define a consistent
    % scaling so the arrow lengths don't vary too much.
    scale_factor = 4e6;

    % Show B field vectors over a regular grid, then overlay colored
    % vectors where there are Qii zero crossings.
    ah = myBox.drawBox( 'show_M', true, 'M_length', 3, 'view', [ 0 90 ] );
    coarse_2D_survey = SurveyField( -10:10, -10:10, 0 );
    myBox.showBfieldVectors( ah, coarse_2D_survey, 'clr', [ .8 .8 .8 ], 'scale_factor', scale_factor );
    Qii_zc_survey = SurveyField( myBox.computeQiiZeros( -10:10 ) );
    myBox.showBfieldVectors( ah, Qii_zc_survey, 'clr', [ .5 .2 .5 ], 'scale_factor', scale_factor );
    text( 0.03, 0.97, '(a)', 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 16 );
    axis( [ -10 10 -10 10 ] );
    exportgraphics( gcf, '../Figures/M_versus_B_directions_Qii_zero_crossings.png', 'ContentType', 'image', 'Resolution', 300 );

    % Show B field vectors over a regular grid, then overlay colored
    % vectors where there are Qij zero crossings (two different kinds).
    ah = myBox.drawBox( 'show_M', true, 'M_length', 3, 'view', [ 0 90 ] );
    coarse_2D_survey = SurveyField( -10:10, -10:10, 0 );
    myBox.showBfieldVectors( ah, coarse_2D_survey, 'clr', [ .8 .8 .8 ], 'scale_factor', scale_factor );
    Qij_zc_survey_x = SurveyField( -10:10, 0, 0 );
    myBox.showBfieldVectors( ah, Qij_zc_survey_x, 'clr', [ .4 .4 .7 ], 'scale_factor', scale_factor );
    Qij_zc_survey_y = SurveyField( 0, -10:10, 0 );
    myBox.showBfieldVectors( ah, Qij_zc_survey_y, 'clr', [ .7 .4 .4 ], 'scale_factor', scale_factor );
    text( 0.03, 0.97, '(b)', 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 16 );
    axis( [ -10 10 -10 10 ] );
    exportgraphics( gcf, '../Figures/M_versus_B_directions_Qij_zero_crossings.png', 'ContentType', 'image', 'Resolution', 300 );

end
if show_3D_Q_and_B

    % Establish a magnetized box and a 3D survey volume.
    myBox = MagBox( [ -2 2 ], [ -1 1 ], [ -.5 .5 ], [ 3 0 -1 ]' );
    survey_volume = SurveyField( linspace(-8,8), linspace(-8,8), linspace(-8,8) );

    % Show Q field.
    fh = BaseTools.tile_figures( myBox.showQfieldContours( 'all', survey_volume, 'title', true ) );
    fh.Position(3:4) = [ 950 1000 ];
    exportgraphics( fh, '../Figures/3D_Q.png', 'ContentType', 'image', 'Resolution', 300 );

    % Show B field.
    fh = BaseTools.tile_figures( myBox.showBfieldContours( 'xyz', survey_volume, 'title', true, 'view', [ 25 15 ] ) );
    exportgraphics( fh, '../Figures/3D_B.png', 'ContentType', 'image', 'Resolution', 300 );

end
if show_Blakely

    % Reproduce finite-cube equivalent of Blakely's figures 4.9 and 4.10.
    fh = MagBox.unit_test('Blakely_Figs_4_9_4_10');
    exportgraphics( fh{1}, '../Figures/Blakely_prism_Fig4-9.png' );
    exportgraphics( fh{2}, '../Figures/Blakely_prism_Fig4-10.png' );
    exportgraphics( fh{3}, '../Figures/Blakely_3D_vert.png' );
    exportgraphics( fh{4}, '../Figures/Blakely_3D_horiz.png' );

end
if show_Bongiolo

    % Reproduce a version of Figure 7 from Bongiolo et al. (2013).
    fh = MagEnsemble.unit_test('Bongiolo7');
    text( fh{1}.CurrentAxes, 0.03, 0.97, '(a)', 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 16 );
    exportgraphics( fh{1}, '../Figures/Bongiolo_Fig7.png' );

    % Show another example with a prism in a 3D rotation.
    fh = MagEnsemble.unit_test('mixed_orientation');
    text( 0.03, 0.97, '(b)', 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 16 );
    exportgraphics( fh, '../Figures/two_prisms_mixed_orientation.png' );

end
if show_wire_plate_limits

    % Show a long wire in 3D.
    L = 100;
    w = 1;
    myBox = MagBox( [ -w/2 w/2 ], [ -L/2 L/2 ], [ -w/2 w/2 ], [ 1 -.2 -.5 ]' );
    survey_volume = SurveyField( linspace(-L,L), linspace(-L,L), linspace(-L,L) );
    fh = BaseTools.tile_figures( myBox.showBfieldContours( 'xyz', survey_volume, 'view', [ -30 30 ], 'title', true ) );
    exportgraphics( fh, '../Figures/3DB_long_wire.png' );

    % Show a wide flat plate in 3D.
    w = 100;
    t = 1;
    myBox = MagBox( [ -w/2 w/2 ], [ -w/2 w/2 ], [ -t/2 t/2 ], [ 1 -.2 -.5 ]' );
    survey_volume = SurveyField( linspace(-w,w), linspace(-w,w), linspace(-w,w) );
    fh = BaseTools.tile_figures( myBox.showBfieldContours( 'xyz', survey_volume, 'view', [ -30 30 ], 'title', true ) );
    exportgraphics( fh, '../Figures/3DB_wide_plate.png' );

end
if show_wire_plate_movie

    % Now make a movie version of this.
    M = [ 3 0 -1 ]';
    keyBoxes = [
        MagBox( [ -1 1 ], [ -1 1 ], [ -0.5 0.5 ], M );
        MagBox( [ -1 1 ], [ -30 30 ], [ -0.5 0.5 ], M );
        MagBox( [ -30 30 ], [ -30 30 ], [ -0.5 0.5 ], M );
        MagBox( [ -30 30 ], [ -1 1 ], [ -0.5 0.5 ], M );
        MagBox( [ -1 1 ], [ -1 1 ], [ -0.5 0.5 ], M );
        ];
    keyBoxes.makeBmovie( '../Figures/3DB_wire_plate_cycle.mp4', 'fixed_moment', 600, 'b', 50, 'duration', 3 );

end
if show_shape_movie

    % Make a movie for a prism whose shape is changing.
    M = [ 3 0 -1 ]';
    keyBoxes = [
        MagBox( [ -1 1 ], [ -1 1 ], [ -0.5 0.5 ], M );
        MagBox( [ -10 10 ], [ -1 1 ], [ -0.5 0.5 ], M );
        MagBox( [ -10 10 ], [ -1 1 ], [ -10 0.5 ], M );
        MagBox( [ -1 1 ], [ -1 1 ], [ -10 0.5 ], M );
        MagBox( [ -1 1 ], [ -1 1 ], [ -0.5 0.5 ], M );
        ];
    keyBoxes.makeBmovie( '../Figures/3DB_shape_cycle.mp4', 'fixed_moment', 200, 'b', 20, 'duration', 3 );
    
end
if show_direction_movie

    % Make a movie for a prism whose magnetization direction is changing.
    wx = 5;
    wy = 2;
    wz = 1;
    keyBoxes = [
        MagBox( [ -wx/2 wx/2 ], [ -wy/2 wy/2 ], [ -wz/2 wz/2 ], [ 3 0 0 ]' );
        MagBox( [ -wx/2 wx/2 ], [ -wy/2 wy/2 ], [ -wz/2 wz/2 ], [ 0 0 3 ]' );
        MagBox( [ -wx/2 wx/2 ], [ -wy/2 wy/2 ], [ -wz/2 wz/2 ], [ -3 0 0 ]' );
        MagBox( [ -wx/2 wx/2 ], [ -wy/2 wy/2 ], [ -wz/2 wz/2 ], [ 0 0 -3 ]' );
        MagBox( [ -wx/2 wx/2 ], [ -wy/2 wy/2 ], [ -wz/2 wz/2 ], [ 0 3 0 ]' );
        MagBox( [ -wx/2 wx/2 ], [ -wy/2 wy/2 ], [ -wz/2 wz/2 ], [ 0 0 3 ]' );
        MagBox( [ -wx/2 wx/2 ], [ -wy/2 wy/2 ], [ -wz/2 wz/2 ], [ 0 -3 0 ]' );
        MagBox( [ -wx/2 wx/2 ], [ -wy/2 wy/2 ], [ -wz/2 wz/2 ], [ 1 -1 -2 ]' );
        MagBox( [ -wx/2 wx/2 ], [ -wy/2 wy/2 ], [ -wz/2 wz/2 ], [ 2 1 -1 ]' );
        MagBox( [ -wx/2 wx/2 ], [ -wy/2 wy/2 ], [ -wz/2 wz/2 ], [ -1 1 2 ]' );
        MagBox( [ -wx/2 wx/2 ], [ -wy/2 wy/2 ], [ -wz/2 wz/2 ], [ -2 0 1 ]' );
        MagBox( [ -wx/2 wx/2 ], [ -wy/2 wy/2 ], [ -wz/2 wz/2 ], [ -1 -1 -1 ]' );
        MagBox( [ -wx/2 wx/2 ], [ -wy/2 wy/2 ], [ -wz/2 wz/2 ], [ 3 0 0 ]' );
        ];
    keyBoxes.makeBmovie( '../Figures/3DB_direction_cycle.mp4', 'fixed_moment', 500, 'b', 40, 'duration', 4 );
    
end
if show_semi_infinite_movies

    % Make a movie showing what happens when one or more of a prism's faces
    % is extended out to infinity.
    M = [ 3 0.5 -1 ]';
    keyBoxes = [
        MagBox( [ -10 0 ], [ -2 2 ], [ -.5 .5 ], M );
        MagBox( [ -100 0 ], [ -2 2 ], [ -.5 .5 ], M );
        MagBox( [ -1000 0 ], [ -2 2 ], [ -.5 .5 ], M );
        MagBox( [ -Inf 0 ], [ -2 2 ], [ -.5 .5 ], M );
        MagBox( [ -Inf 10 ], [ -2 2 ], [ -.5 .5 ], M );
        MagBox( [ -Inf 100 ], [ -2 2 ], [ -.5 .5 ], M );
        MagBox( [ -Inf 1000 ], [ -2 2 ], [ -.5 .5 ], M );
        MagBox( [ -Inf Inf ], [ -2 2 ], [ -.5 .5 ], M );
        MagBox( [ -1000 Inf ], [ -2 2 ], [ -.5 .5 ], M );
        MagBox( [ -100 Inf ], [ -2 2 ], [ -.5 .5 ], M );
        MagBox( [ -10 Inf ], [ -2 2 ], [ -.5 .5 ], M );
        MagBox( [ -10 100 ], [ -2 2 ], [ -.5 .5 ], M );
        MagBox( [ -10 10 ], [ -2 2 ], [ -.5 .5 ], M );
        MagBox( [ -10 0 ], [ -2 2 ], [ -.5 .5 ], M );
        ];
    keyBoxes.makeBmovie( '../Figures/semi-infinite-cycle-both-sides.mp4', 'b', 50, 'duration', 3 );

end

