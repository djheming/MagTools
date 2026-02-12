% Reproduce figures for the paper accompanying this software (Hemingway,
% RASTI, 2026).
%
% SCRIPT USAGE: This script is divided into sections. Instead of running
% the whole script at once, the recommended approach is to click into a
% section and press Cmd+Enter (Mac) or Ctrl+Enter (PC) to run one section
% at a time. Each section corresponds to one figure in the paper (or a few
% if they're closely related).
%
% If you want figures to be saved, specify a valid destination folder below
% (modify the line that sets figs_folder = []).
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


%% Step 0: Initialize.

% Before doing anything, make sure MagTools and BaseTools (a library for
% lower level functions) are on your Matlab path. You can do this by
% running the setup.m script.  
setup;

% Clean up previous figures, if desired.
% close all; % Uncomment this line if you want to clear previous figures each time you run this script.

% If you want the figures saved, specify a destination folder here. 
% If this is empty, figures are still displayed but not saved.
figs_folder = [];


%% Figures 1 and 6: 3D Box Geometry
% Illustrates the important vector definitions for the general case of a
% rectangular prism. Also illustrates the lines and solid angles (pyramids)
% discussed in the geometric interpretation of the Qii and Qij functions.

% Establish a magnetized box and define an example evaluation point.
myBox = MagBox( [ 1.5 3.5 ], [ 0.5 1.5 ], [ 0.3 1.3 ], [ 1 0 0 ]' );
p = [ 2.45 -1 0.9 ]';

% Show the box in 3D, illustrating the evaluation point.
ah = myBox.drawBox( 'p', p, 'axes', true, 'vertices', true, 'show_M', false, 'show_p_vector', true, 'label_p', true, 'verts_from_origin', { '122', '211' }, 'view', [ 39 47 ] );
fh = ah.Parent;

% Save orthographic and oblique views.
if exist( 'figs_folder', 'var' ) && ~isempty( figs_folder )
    view( ah, [ 0 90 ] ); % xy
    exportgraphics( fh, [ figs_folder '/3D_box_geometry_xy.png' ], 'ContentType', 'image', 'Resolution', 600 );
    view( ah, [ 0 0 ] ); % xz
    exportgraphics( fh, [ figs_folder '/3D_box_geometry_xz.png' ], 'ContentType', 'image', 'Resolution', 600 );
    view( ah, [ 90 0 ] ); % yz
    exportgraphics( fh, [ figs_folder '/3D_box_geometry_yz.png' ], 'ContentType', 'image', 'Resolution', 600 );
    view( ah, [ 39 47 ] );
    exportgraphics( fh, [ figs_folder '/Figure1_3D_box_geometry_oblique.png' ], 'ContentType', 'image', 'Resolution', 600 );
end

% Make a box and define a couple of different evaluation points.
myBox = MagBox( [ 1.5 3.5 ], [ 0.5 1.5 ], [ -1 0 ], [ 1 0 0 ]' );
p = [
    2.45 -1 0.5 % Chosen such that Qii<0 and Qij is nearly zero
    4.2 -1 -0.4 % Chosen such that Qii is nearly zero and Qij is non-zero
    ]';
num_points = size(p,2);

% For each evaluation point, show the box with shaded solid angles on
% the +i and -i faces.
for k = 1 : num_points
    ah = myBox.drawBox( 'p', p(:,k), 'vertices', true, 'show_M', false, 'show_p_vector', true, 'label_p', true, 'shading', true, 'axes', true, 'axislabels', 'ijk', 'rlabels', false, 'alphas_and_fs', true, 'view', [ 20 60 ] );
    fh = ah.Parent;
    axis( [ 0 4.5 -1 2 -1.5 2 ] );
    if exist( 'figs_folder', 'var' ) && ~isempty( figs_folder )
        exportgraphics( fh, [ figs_folder '/Figure6_box_geometry_with_lines_and_angles_' num2str(k) '.png' ], 'Resolution', 600 );
    end
end


%% Figure 2: 3D Prism Coordinate Transformation.
% Illustrates the relationship between the source-aligned coordinate system (S) 
% and a more general analysis coordinate system (A). 

% Define a prism that is a tad askew.
roll = 10;
pitch = -20;
yaw = 30;
R_AS = BaseTools.rpy2rot( roll, pitch, yaw );
v_AS = [ 2 -.5 0 ]';
M = [ 0 1 0 ]'; % This is in the source coordinate frame.
myBox = MagBox( [ 0 2 ], [ 0 1 ], [ 0 0.6 ], M, v_AS, R_AS, Mframe='S' );

% Show box with evaluation point.
p = [ 2.5 1.5 0 ]';
ah = myBox.drawBox( 'axes', true, 'p', p, 'show_M', false, 'label_origin', true, 'show_p_vector', true, 'label_p_vector', true, 'view', [ 0 90 ] );
axis( ah, 'off' );
% camlight;
if exist( 'figs_folder', 'var' ) && ~isempty( figs_folder )
    exportgraphics( gcf, [ figs_folder '/Figure2_3D_coordinate_change.png' ], 'ContentType', 'image', 'Resolution', 600 );
end


%% Figure 3: 2D Box Geometry
% Illustrates geometric interpretation of Qii and Qij functions for the 2D
% case.

% Create a 2D prism (infinitely long in one dimension).
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
    ah = myBox.drawBox( 'p', p(:,k), 'axislabels', 'ijk', 'show_M', false, 'label_p', true, 'vertices', true, 'alphas_and_fs', true, 'zerocrossings', true, 'view', [ 0 90 ] );
    if exist( 'figs_folder', 'var' ) && ~isempty( figs_folder )
        exportgraphics( ancestor(ah,'figure'), sprintf('%s/Figure3_2D_box_geometry_p(%s,%s).png', figs_folder, num2str(p(1,k)), num2str(p(2,k)) ), 'Resolution', 600 );
    end
end


%% Figures 4 and 5: Q and B Fields for 2D Prism
% Illustration of the four components of Q for a 2D prism.
% Additionally shows B-field vectors on a 2D cross-section.
% The field vectors along the zero-crossings are highlighted.

% For the 2D problem, make a box that has infinite extent in the
% z-direction.
myBox = MagBox( [ -3 3 ], [ -1 1 ], [ -Inf Inf ], [ 3 1 0 ]' );
fine_2D_survey = SurveyField( linspace(-12,12), linspace(-12,12), 0 );

% Show Q field.
fh = BaseTools.tileFigures( myBox.showQfieldContours( { 'xx', 'xy'; 'yx' 'yy' }, fine_2D_survey, 'title', true ) );
fh.Position(4) = 920;
if exist( 'figs_folder', 'var' ) && ~isempty( figs_folder )
    exportgraphics( fh, [ figs_folder '/Figure4_2D_Q.png' ], 'ContentType', 'image', 'Resolution', 600 );
end

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
fha = gcf;

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
fhb = gcf;

% Combine the two figure panels.
fh = BaseTools.tileFigures( [ fha fhb ] );
if exist( 'figs_folder', 'var' ) && ~isempty( figs_folder )
    exportgraphics( fh, [ figs_folder '/Figure5_M_versus_B_directions_zero_crossings.png' ], 'ContentType', 'image', 'Resolution', 600 );
end


%% Figures 7 and 8: 3D Q and B field contours.
% For a finite 3D prism, shows all nine components of the Q field and the
% three components of a B-field for an assumed magnetization direction.

% Assume some magnetization direction.
M = [ 3 0 -1 ]';

% Establish a magnetized box and a 3D survey volume.
myBox = MagBox( [ -2 2 ], [ -1 1 ], [ -.5 .5 ], M );
survey_volume = SurveyField( linspace(-8,8), linspace(-8,8), linspace(-8,8) );

% Show Q field.
fh = BaseTools.tileFigures( myBox.showQfieldContours( 'all', survey_volume, 'title', true ) );
fh.Position(3:4) = [ 950 1000 ];
if exist( 'figs_folder', 'var' ) && ~isempty( figs_folder )
    exportgraphics( fh, [ figs_folder '/Figure7_3D_Q.png' ], 'ContentType', 'image', 'Resolution', 600 );
end

% Show B field.
fh = BaseTools.tileFigures( myBox.showBfieldContours( 'xyz', survey_volume, 'title', true, 'view', [ 25 15 ] ) );
if exist( 'figs_folder', 'var' ) && ~isempty( figs_folder )
    exportgraphics( fh, [ figs_folder '/Figure8_3D_B.png' ], 'ContentType', 'image', 'Resolution', 600 );
end


%% Figures A1-A3: Reproduce Blakely (1995) Figures 4.9 and 4.10.

fh = MagBox.unit_test('Blakely_Figs_4_9_4_10');
if exist( 'figs_folder', 'var' ) && ~isempty( figs_folder )
    exportgraphics( fh{1}, [ figs_folder '/FigureA1_Blakely_prism_Fig4-9.png' ] );
    exportgraphics( fh{2}, [ figs_folder '/FigureA2_Blakely_prism_Fig4-10.png' ] );
    exportgraphics( fh{3}, [ figs_folder '/FigureA3_Blakely_3D_vert.png' ] );
    exportgraphics( fh{4}, [ figs_folder '/FigureA3_Blakely_3D_horiz.png' ] );
end


%% Figures A4 and A5: 3D Visualizations of Long Wire and Flat Plate Cases

% Show a long wire in 3D.
L = 100;
w = 1;
myBox = MagBox( [ -w/2 w/2 ], [ -L/2 L/2 ], [ -w/2 w/2 ], [ 1 -.2 -.5 ]' );
survey_volume = SurveyField( linspace(-L,L), linspace(-L,L), linspace(-L,L) );
fh = BaseTools.tileFigures( myBox.showBfieldContours( 'xyz', survey_volume, 'view', [ -30 30 ], 'title', true ) );
if exist( 'figs_folder', 'var' ) && ~isempty( figs_folder )
    exportgraphics( fh, [ figs_folder '/FigureA4_3DB_long_wire.png' ] );
end

% Show a wide flat plate in 3D.
w = 100;
t = 1;
myBox = MagBox( [ -w/2 w/2 ], [ -w/2 w/2 ], [ -t/2 t/2 ], [ 1 -.2 -.5 ]' );
survey_volume = SurveyField( linspace(-w,w), linspace(-w,w), linspace(-w,w) );
fh = BaseTools.tileFigures( myBox.showBfieldContours( 'xyz', survey_volume, 'view', [ -30 30 ], 'title', true ) );
if exist( 'figs_folder', 'var' ) && ~isempty( figs_folder )
    exportgraphics( fh, [ figs_folder '/FigureA5_3DB_wide_plate.png' ] );
end


%% Figure A6: Bongiolo replication and arbitrary 3D rotation example

% Reproduce a version of Figure 7 from Bongiolo et al. (2013).
fh = MagEnsemble.unit_test('Bongiolo7');
text( fh{1}.CurrentAxes, 0.03, 0.97, '(a)', 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 16 );
if exist( 'figs_folder', 'var' ) && ~isempty( figs_folder )
    exportgraphics( fh{1}, [ figs_folder '/FigureA6a_Bongiolo_Fig7.png' ] );
end

% Show another example with a prism in a 3D rotation.
fh = MagEnsemble.unit_test('mixed_orientation');
text( 0.03, 0.97, '(b)', 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 16 );
if exist( 'figs_folder', 'var' ) && ~isempty( figs_folder )
    exportgraphics( fh, [ figs_folder '/FigureA6b_two_prisms_mixed_orientation.png' ] );
end


%% Supplementary Video: Direction Movie
% Movie for a prism whose magnetization direction is changing.

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
if exist( 'figs_folder', 'var' ) && ~isempty( figs_folder )
    keyBoxes.makeBmovie( [ figs_folder '/Video1_3DB_direction_cycle.mp4' ], 'fixed_moment', 500, 'b', 40, 'duration', 4 );
else
    keyBoxes.previewBmovie( 'fixed_moment', 500, 'b', 40 );
end


%% Supplementary Video: Shape Movie
% Movie for a prism whose shape is changing.

M = [ 3 0 -1 ]';
keyBoxes = [
    MagBox( [ -1 1 ], [ -1 1 ], [ -0.5 0.5 ], M );
    MagBox( [ -10 10 ], [ -1 1 ], [ -0.5 0.5 ], M );
    MagBox( [ -10 10 ], [ -1 1 ], [ -10 0.5 ], M );
    MagBox( [ -1 1 ], [ -1 1 ], [ -10 0.5 ], M );
    MagBox( [ -1 1 ], [ -1 1 ], [ -0.5 0.5 ], M );
    ];
if exist( 'figs_folder', 'var' ) && ~isempty( figs_folder )
    keyBoxes.makeBmovie( [ figs_folder '/Video2_3DB_shape_cycle.mp4' ], 'fixed_moment', 200, 'b', 20, 'duration', 3 );
else
    keyBoxes.previewBmovie( 'fixed_moment', 200, 'b', 20 );
end


%% Supplementary Video: Wire Plate Movie
% Movie for a prism that transitions between a long wire and a wide plate.

M = [ 3 0 -1 ]';
keyBoxes = [
    MagBox( [ -1 1 ], [ -1 1 ], [ -0.5 0.5 ], M );
    MagBox( [ -1 1 ], [ -30 30 ], [ -0.5 0.5 ], M );
    MagBox( [ -30 30 ], [ -30 30 ], [ -0.5 0.5 ], M );
    MagBox( [ -30 30 ], [ -1 1 ], [ -0.5 0.5 ], M );
    MagBox( [ -1 1 ], [ -1 1 ], [ -0.5 0.5 ], M );
    ];
if exist( 'figs_folder', 'var' ) && ~isempty( figs_folder )
    keyBoxes.makeBmovie( [ figs_folder '/Video3_3DB_wire_plate_cycle.mp4' ], 'fixed_moment', 600, 'b', 50, 'duration', 3 );
else
    keyBoxes.previewBmovie( 'fixed_moment', 600, 'b', 50 );
end


%% Supplementary Video: Testing the Semi-infinite Cases
% Movie of a prism that is extended out to infinity on one side, then both
% sides, and then returns to having finite length.

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
if exist( 'figs_folder', 'var' ) && ~isempty( figs_folder )
    keyBoxes.makeBmovie( [ figs_folder '/Video4_semi-infinite-cycle-both-sides.mp4' ], 'b', 50, 'duration', 3 );
else
    keyBoxes.previewBmovie( 'b', 50 );
end


