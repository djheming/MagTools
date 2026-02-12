function B = computeBfield( thisMagBox, p, varargin )

% B = computeBfield( thisMagBox, p ) compute magnetic field (B),
% evaluated at an array of positions (p), arising from the supplied
% magnetized rectangular prism (thisMagBox). p is a 3xN matrix whose rows
% represent the x, y, z, components, and whose columns represent N
% different positions. 
% 
% These points are always expressed in the "analysis" coordinate system,
% which may or may not differ from the "source" coordinate system. If they
% differ, some internal logic is required to carry out the necessary
% transformations within the computeQfield function. 

% The output is a 3xN matrix whose rows represent the Bx, By, Bz components
% and whose columns represent the N different vector field measurements.
% Again, this is always in the "analysis" coordinate system.
%
%
%   Disclaimer: This code is provided as-is, has been tested only very
%   informally, and may not always behave as intended. I find it useful for
%   my own work, and I hope you will too, but I make no guarantees as to
%   the accuracy or robustness of this code. This code is also actively
%   under development and future versions may not be backward compatible.
%
%   Created: 2018-06-07
%   Doug Hemingway
%   University of California, Berkeley
%
%   Updated: 2025-04-14
%   Doug Hemingway (douglas.hemingway@utexas.edu)
%   University of Texas at Austin
%


% Parse input arguments.
args = BaseTools.argarray2struct( varargin );

% Make sure the input position grid looks right.
if size( p, 1 ) ~= 3
    error( 'Input position vector array requires 3 rows' );
end
Np = size( p, 2 );

% Obtain the Q matrix for each of the observation points (p). The result
% here will be a 3x3xN matrix expressing the Q transformation in the
% analysis coordinate system (A). We can skip this if the user has supplied a
% Q with the right dimensions.
if isfield( args, 'Q' ) && ~isempty( args.Q ) && all( size(args.Q)==[ 3 3 Np ] )
    Q = args.Q;
else
    Q = thisMagBox.computeQfield( p );
end

% Combine this with the magnetization vector to get the B field. 
B = MagSource.mu_naught_over_4pi * pagemtimes( Q, thisMagBox.M );

% Squeeze out the singleton dimension to leave only a 3xN matrix.
B = squeeze(B);


