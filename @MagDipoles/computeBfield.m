function B = computeBfield( theseDipoles, p, varargin )

% B = computeBfield( theseDipoles, p ) compute the magnetic flux density
% (B) evaluated at an array of positions (p). For each position, B is a 3x1
% vector representing the vector B field.
% 
% The input array p is a 3xN matrix whose rows represent the x, y, z,
% components, and whose columns represent N different positions. 
% The output is a 3xN array. 
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
%   Updated: 2025-06-06
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

% Before carrying out the B = mu_naught_over_4pi * Q * M, we must make sure
% that each dipole has some effective volume (otherwise, M is not defined).
% If no effective volume (dV) is yet set, we can just set them all to have
% unit volume.
if isempty( theseDipoles.dV )
    theseDipoles.dV = 1;
end
   
% First, obtain the full Q, allowing for the possibility of distinct
% magnetizations for each dipole. This will be a 3 x 3 x Np x Nq array. We
% can skip this if the user has supplied a Q with the right dimensions.
if isfield( args, 'Q' ) && ~isempty( args.Q ) && all( size(args.Q)==[ 3 3 Np theseDipoles.Nq ] )
    Q = args.Q;
else
    Q = theseDipoles.computeQfield( p, 'merge_dipoles', false );
end

% We want to multiply this by M, again allowing for the possibility of distinct
% magnetizations for each dipole. To do this, we need to reshape M to have
% shape 3 x 1 x 1 x Nq.
M_reshaped = reshape( theseDipoles.M, 3, 1, 1, [] );

% Now matrix multiply each of the Nq dipoles at each of the Np points. This
% gives us a 3 x 1 x Np x Nq array.
V = pagemtimes( Q, M_reshaped );

% Then we sum across the Nq dipoles (the 4th dimension). At the same time,
% we can multiply by mu_naught_over_4pi.
B = MagSource.mu_naught_over_4pi * sum( V, 4 );

% Finally, we need to reshape to get the desired 3 x 1 x Np array.
B = reshape( B, 3, [] );


