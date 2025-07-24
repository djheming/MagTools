function Q = computeQfield( theseDipoles, p, varargin )

% Q = computeQfield( theseDipoles, p ) compute the Q matrix 
% evaluated at an array of positions (p). For each position, Q is a 3x3
% matrix that depends on the location and effective volume of each dipole.
% If the dipoles do not have associated effective volumes (dV), this
% won't work. 
% 
% The input array p is a 3xN matrix whose rows represent the x, y, z,
% components, and whose columns represent N different positions. 
% The output is a 3x3xN matrix. The matrices can be used to relate the
% magnetization vector (M) to the resulting B-field. 
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


% Make sure the dipoles have associated effective volumes.
if isempty(theseDipoles.dV) || any( theseDipoles.dV==0 )
    error( 'dV must be specified before we can compute Q.' );
end

% Make sure the input position grid looks right.
if size( p, 1 ) ~= 3
    error( 'Input position vector array requires 3 rows.' );
end

% Parse additional arguments.
args = BaseTools.argarray2struct( [ { 'merge_dipoles', true 'qtol', theseDipoles.dV/5 } varargin ] );

% Compute r vectors. We are transposing the q vectors such that the
% resulting r vector components have size Nq x N where Nq is the number of
% dipoles in the array.
rx = p(1,:) - theseDipoles.q(1,:)';
ry = p(2,:) - theseDipoles.q(2,:)';
rz = p(3,:) - theseDipoles.q(3,:)';

% Compute six independent elements of the Q matrix. Each of these has dimensionality 1x1xN.
Qxx = compute_Qii( rx, ry, rz, theseDipoles.dV' );
Qxy = compute_Qij( rx, ry, rz, theseDipoles.dV' );
Qxz = compute_Qij( rx, rz, ry, theseDipoles.dV' );
Qyx = Qxy;
Qyy = compute_Qii( ry, rx, rz, theseDipoles.dV' );
Qyz = compute_Qij( ry, rz, rx, theseDipoles.dV' );
Qzx = Qxz;
Qzy = Qyz;
Qzz = compute_Qii( rz, ry, rx, theseDipoles.dV' );

% Package up the results in a 3 x 3 x Np x Nq matrix.
Q = [
    Qxx Qxy Qxz
    Qyx Qyy Qyz
    Qzx Qzy Qzz
    ];

% In some cases, the user may wish to merge the contributions of the
% separate dipoles (e.g., when plotting the Q fields alone).
if args.merge_dipoles
    Q = sum(Q,4);
end

% Things get weird if you're too close to the dipoles, so we need to NaN
% out a small volume around each dipole. By default, we'll define the
% exclusion zone with some radius qtol, which is proportional to dV.
d = pdist2( p', theseDipoles.q' ); % returns Np × Nq matrix of distances
tooclose = any( d < args.qtol, 2 ); % logical index (true if too close to dipoles)
Q(:,:,tooclose) = NaN;



function Qii = compute_Qii( ri, rj, rk, dV )

% The input r vector elements are each Nq x Np arrays.
% dV may be a scalar or an Nq x 1 array.
rsq = ( ri.^2 + rj.^2 + rk.^2 );
Qii = ( ( 2*ri.^2 - rj.^2 - rk.^2 )./rsq.^(5/2) ).*dV;

% At this point, Qii is a Nq x Np array. We want the output to be a 1 x 1 x
% Np x Nq array so we'll have to permute the dimensions.
Qii = permute( Qii, [ 3 4 2 1 ] );




function Qij = compute_Qij( ri, rj, rk, dV )

% The input r vector elements are each Nq x N arrays.
% dV may be a scalar or an Nqx1 array.
rsq = ( ri.^2 + rj.^2 + rk.^2 );
Qij = ( ( 3 * ri .* rj )./rsq.^(5/2) ).*dV;

% At this point, Qij is a Nq x Np array. We want the output to be a 1 x 1 x
% Np x Nq array so we'll have to permute the dimensions.
Qij = permute( Qij, [ 3 4 2 1 ] );





