function B = computeBfield( thisMagBox, p, varargin )

% B = computeBfield( thisMagBox, p ) compute magnetic field (B),
% evaluated at an array of positions (p), arising from the supplied
% magnetized rectangular prism (thisMagBox). p is a 3xN matrix whose rows
% represent the x, y, z, components, and whose columns represent N
% different positions. The output is likewise a 3xN matrix whose rows
% represent the Bx, By, Bz components and whose columns represent the N
% different vector field measurements.
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

% First, we have to make sure to express the evaluation points (p) in the
% magnetized prism's own preferred "source coordinate system" (S).
pS = thisMagBox.R_SA * ( p - thisMagBox.v_AS );

% Obtain the Q matrix for each of the observation points (p). The result
% here will be a 3x3xN matrix expressing the Q transformation in the
% source coordinate system (S). We can skip this if the user has supplied a
% Q with the right dimensions.
if isfield( args, 'Q' ) && ~isempty( args.Q ) && all( size(args.Q)==[ 3 3 Np ] )
    Q_S = args.Q;
else
    Q_S = thisMagBox.computeQfield( pS );
end

% Combine this with the magnetization vector to get the B field. The way we
% do this depends on whether the magnetization vector is expressed in the
% analysis coordinate system or the source coordinate system.
switch thisMagBox.Mframe
    case 'S'
        Q_A = pagemtimes( thisMagBox.R_AS, Q_S );
    case 'A'
        Q_A = pagemtimes( thisMagBox.R_AS, pagemtimes( Q_S, thisMagBox.R_SA ) );
    otherwise
        error( 'Unrecognized coordinate frame for magnetization vector' );
end
B = MagSource.mu_naught_over_4pi * pagemtimes( Q_A, thisMagBox.M );

% Squeeze out the singleton dimension to leave only a 3xN matrix.
B = squeeze(B);


