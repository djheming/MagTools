function Q = computeQfield( thisMagBox, p, varargin )

% Q = computeQfield( thisMagBox, p ) compute the Q matrix 
% evaluated at an array of positions (p). 
% 
% For each position, Q is a 3x3 matrix that depends on the geometry of the
% magnetized box (thisMagBox) and the position, p. 
% 
% The input array p is a 3xN matrix whose rows represent the x, y, z,
% components, and whose columns represent N different positions. 
%
% The output is a 3x3xN matrix. The matrices can be used to relate the
% magnetization vector (M) to the resulting B-field (see computeBfield
% for MagBox). 
%
% IMPORTANT: Internally, the terms in the Q matrix are calculated from the
% r vectors, which represent the evaluation points relative to the prism's
% eight vertices. To work properly, the evaluation points MUST be expressed
% in the prism's own coordinate frame (the "source" frame). By default,
% however, the external logic wants the Q matrix to be expressed in the
% "analysis" coordinate frame, so some internal manipulations are required
% at the front and back end of this function. The user has the option of
% overriding those bookend transformations (by specifying Qframe='S') in
% order to get direct access to the Q matrix in its raw "source" frame
% form.
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
% By default, assume the calling function wants to use the "analysis" frame.
args = BaseTools.argarray2struct( varargin, { 'Qframe', 'A' } ); 
if strcmp(args.Qframe,'A')
    % Here, the input evaluation points (p) are expressed in the "analysis"
    % frame. We have to transform them into the "source" frame before
    % we can compute Q.
    pS = thisMagBox.R_SA * ( p - thisMagBox.v_AS );
elseif strcmp(args.Qframe,'S')
    % User is already working in the "source" frame.
    pS = p;
else
    error( 'Unrecognized coordinate frame (%s)', args.Qframe );
end

% Make sure the input position grid looks right.
if size( pS, 1 ) ~= 3
    error( 'Input position vector array requires 3 rows' );
end

% Compute integration limit vectors.
rx1 = pS(1,:) - min(thisMagBox.vx);
rx2 = pS(1,:) - max(thisMagBox.vx);
ry1 = pS(2,:) - min(thisMagBox.vy);
ry2 = pS(2,:) - max(thisMagBox.vy);
rz1 = pS(3,:) - min(thisMagBox.vz);
rz2 = pS(3,:) - max(thisMagBox.vz);

% Compute six independent elements of the Q matrix. Each of these has dimensionality 1x1xN.
Qxx = MagBox.computeQii( [ rx1; rx2 ], [ ry1; ry2 ], [ rz1; rz2 ] );
Qxy = MagBox.computeQij( [ rx1; rx2 ], [ ry1; ry2 ], [ rz1; rz2 ] );
Qxz = MagBox.computeQij( [ rx1; rx2 ], [ rz1; rz2 ], [ ry1; ry2 ] );
Qyx = Qxy;
Qyy = MagBox.computeQii( [ ry1; ry2 ], [ rx1; rx2 ], [ rz1; rz2 ] );
Qyz = MagBox.computeQij( [ ry1; ry2 ], [ rz1; rz2 ], [ rx1; rx2 ] );
Qzx = Qxz;
Qzy = Qyz;
Qzz = MagBox.computeQii( [ rz1; rz2 ], [ ry1; ry2 ], [ rx1; rx2 ] );

% Package up the results in a 3x3xN matrix.
Q = [
    Qxx Qxy Qxz
    Qyx Qyy Qyz
    Qzx Qzy Qzz
    ];

% These formulas don't work inside the box so we need to NaN out those parts.
insidebox = pS(1,:)>=min(thisMagBox.vx) & pS(1,:)<=max(thisMagBox.vx) & ...
    pS(2,:)>=min(thisMagBox.vy) & pS(2,:)<=max(thisMagBox.vy) & ...
    pS(3,:)>=min(thisMagBox.vz) & pS(3,:)<=max(thisMagBox.vz);
Q(:,:,insidebox) = NaN;

% If the calling function is working in the "analysis" frame, we have to
% rotate the resulting Q from the (S) frame back to the (A) frame. If not,
% then by process of elimination (i.e., we already checked the
% possibilities above), the calling function wants Q in the "source" frame,
% which it already is, so we do nothing.
if strcmp(args.Qframe,'A')
    Q = pagemtimes( thisMagBox.R_AS, pagemtimes( Q, thisMagBox.R_SA ) );
end
