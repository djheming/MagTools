% This script is intended to serve as a set of basic examples of how to use
% MagTools.
%
% MagTools is a set of tools for modeling magnetized source bodies in
% Matlab. The software is organized in an object-oriented framework with
% several interacting classes.  
%
% Typical usage involves defining a magnetized object and a region over
% which to display the resulting field and then using functions like
% showBfieldContours or showBfieldVectors to visualize the resulting
% magnetic field. See below for details.
%
% SCRIPT USAGE: This script is divided into sections. Instead of running
% the whole script at once, the recommended approach is to click into a
% section and press Cmd+Enter (Mac) or Ctrl+Enter (PC) to run one example
% at a time.
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
%   2026-02-02
%


%% Step 0: Initialize.

% Before doing anything, make sure MagTools and BaseTools (a library for
% lower level functions) are on your Matlab path. You can do this by
% running the setup.m script.  
setup;

% Clean up previous figures, if desired.
% close all; % Uncomment this line if you want to clear previous figures each time you run this script.



%% Example 1: A simple dipole.
% This section introduces the use of the MagDipoles class.

% Example 1.1: Define your first magnetized object. Let's make it a simple
% dipole. We'll need to define both a location for the dipole in Cartesian
% space (q) and a magnetic moment (m). Both of these are 3x1 vectors.
% Once we define those, we can pass them into the constructor for the
% MagDipoles class, which returns our object. 
q = [ 3; 2; -1 ]; % Could be any distance units you like
m = [ 4; 0; -1 ]; % This is assumed to be in Am^2.
myDipole = MagDipoles( q, m ) % Dropping the semicolon to reveal the contents of this object.

% Example 1.2: Visualize the magnetic field around this dipole. Let's start
% with just the x-component. The function showBfieldContours is a method
% that works for any type of magnetized object in MagTools. This could be a
% dipole, an array of dipoles, a prism, or an ensemble of multiple types of
% magnetized objects. For now, we'll just start with our single dipole.
myDipole.showBfieldContours('x');

% Example 1.3: Make a more explicit choice about how you want to visualize
% the B-field. In Example 1.2, you did not explicitly state where you want
% to visualize the field so the showBfieldContours function internally
% decided to focus on a region that is 2x2x2 distance units in volume,
% centered on your dipole. Here, we will instead define a specific 1D
% survey transect using the SurveyField class. We'll make it slightly above
% the z=0 plane. In its most typical usage, the SurveyField constructor
% accepts three arguments, defining the x, y, and z grid coordinates,
% respectively. 
xvals = linspace(-10,10,101); % gives a finely spaced set of points along the x-axis
yvals = 0;
zvals = 0;
mySurveyTransect = SurveyField( xvals, yvals, zvals );
myDipole.showBfieldContours('x',mySurveyTransect); 
% Here, we have used a SurveyField object to explicitly tell
% showBfieldContours the region we're interested in.

% Example 1.4: If you ask to view multiple components (e.g., 'xyz') along a
% 1D transect, you get a multi-line plot.
myDipole.showBfieldContours('xyz',mySurveyTransect); 

% Example 1.5: Let's do the same again but now with a 2D plane.
mySurveyPlane = SurveyField( linspace(-10,10,101), linspace(-10,10,101), 0 );
myDipole.showBfieldContours('x',mySurveyPlane);
% Notice that you can freely rotate the figure to see the field from
% different angles. Notice also that, in addition to the high resolution
% shading, you'll see an illustration of vectors, but over a somewhat
% coarser grid. This coarse grid is determined automatically.

% Example 1.6: If you want to see the vectors without the shading, or you
% want more explicit control of where you see the vectors, you'll first
% have to define a coarser grid, and then use the showBfieldVectors
% function. 
coarseSurveyPlane = SurveyField( -10:0.5:10, -10:0.5:10, 0 );
myDipole.showBfieldVectors(coarseSurveyPlane); % No need to specify the component(s) since you'll be seeing the whole 3D vector at each point.

% Example 1.7: For the contour plots, we can also ask for multiple
% components at the same time. The tileFigures function allows us to
% package them all into one figure. 
BaseTools.tileFigures( myDipole.showBfieldContours('xyz',mySurveyPlane) );

% Example 1.8: Let's do the same again but for a 3D volume.
mySurveyVolume = SurveyField( linspace(-5,10), linspace(-5,5), linspace(0,4) );
BaseTools.tileFigures( myDipole.showBfieldContours('xyz',mySurveyVolume) );



%% Example 2: Rectangular prism.
% This section introduces the use of the MagBox class.
% This example builds on the previous example, so execute Example 1 first.

% Example 2.1: Let's move on to a magnetized prism. We'll need to define
% its extent (i.e., min and max) in each of the x, y, and z directions.
% We'll also need a magnetization vector, M.
M = [ 4; 0; -1 ]; % This is assumed to be in A/m.
xrng = [ -2 4 ];
yrng = [ 1 3 ];
zrng = [ -2 -5 ]; % The order is not important.
myBox = MagBox( xrng, yrng, zrng, M ) % Dropping the semicolon to reveal the contents of this object.
BaseTools.tileFigures( myBox.showBfieldContours( 'xyz' ) );

% Example 2.2: Inspect the field on our 2D plane. 
BaseTools.tileFigures( myBox.showBfieldContours( 'xyz', mySurveyPlane ) );

% Example 2.3: Maybe you want to see the vectors without the shading?
% Again, you're going to want a coarser grid for this.
myBox.showBfieldVectors( coarseSurveyPlane );



%% Example 3: Multiple prisms.
% This section introduces the MagEnsemble class.

% Example 3.0: Assemble a collection of prisms with different positions and
% orientations. Here, we'll make two identical boxes but the second will be
% translated and rotated with respect to the first. We will always be
% viewing the resulting field structure in terms of a common "analysis"
% coordinate frame. 

% Box 1: Standard definition, expressed in the "analysis" coordinate frame.
M = [ 3 0 0 ]'; 
xrng = [ -2 4 ];
yrng = [ 1 3 ];
zrng = [ -2 -5 ]; % The order is not important.
box1 = MagBox( xrng, yrng, zrng, M );

% Box 2: Use a vector and rotation matrix to define how this prism's "source"
% coordinate frame relates to the "analysis" coordinate frame.
roll = 5; pitch = -10; yaw = 40;
R_AS = BaseTools.rpy2rot( roll, pitch, yaw ); % Rotation built from roll, pitch, yaw angles (degrees).
v_AS = [ 10; -10; 0 ]; % Vector offset for origin of second prism.
box2 = MagBox( xrng, yrng, zrng, M, v_AS, R_AS ); % Here, M is still expressed in terms of the "analysis" frame.

% Combine the two boxes into an ensemble.
myEnsemble = MagEnsemble( [ box1; box2 ] ) % Dropping the semicolon to reveal the contents of this object.

% Example 3.1: Examine the resulting field on a wide survey plane.
wideSurvey = SurveyField( linspace(-20,20,101), linspace(-20,20,101), 0 );
BaseTools.tileFigures( myEnsemble.showBfieldContours( 'xyz', wideSurvey ) );

% Example 3.2: Here we'll define a third box where the magnetization is
% defined not in terms of the "analysis" coordinate frame but rather in a
% frame that maintains its orientation with the prism (i.e., in the
% "source" frame). 
box3 = MagBox( xrng, yrng, zrng, M, v_AS, R_AS, Mframe='S' );
myEnsemble = MagEnsemble( [ box1; box3 ] )
BaseTools.tileFigures( myEnsemble.showBfieldContours( 'xyz', wideSurvey ) );

% Example 3.3: Let's thicken the survey into a volume.
wideSurvey.zv = 0:.1:5;
BaseTools.tileFigures( myEnsemble.showBfieldContours( 'xyz', wideSurvey ) );



%% Looking for more?

% To reproduce the figures found in the accompanying publication
% (Hemingway, RASTI, 2026), run reproduce_paper_figures.m. For details see
% the contents of that file.

% For many more examples, including more advanced testing and validation,
% refer to the unit_test functions located at the bottom of each of the class
% definition files (e.g., MagBox.m, MagDipoles.m, MagEnsemble.m).

