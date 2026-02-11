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
% IMPORTANT: by default, the positions defined by p must be expressed
% in the magnetized prism's own "source" coordinate system.
% Transformations between an external "analysis" coordinate system and
% this "source" coordinate system are generally handled externally.
% However, if the user specifies Qframe='A', the transformation is handled
% internally. The output is then Q^{(A)}. To get B^{(A)}, this must be
% multipled with M^{(A)}. 
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
args = BaseTools.argarray2struct( varargin, { 'Qframe', 'S' } ); 
if strcmp(args.Qframe,'A')
    % By default, we work in the "source" coordinate frame. But in this case,
    % the user wants us to take inputs and return outputs in the "analysis"
    % coordinate frame instead. First, we transform the points into the
    % "source" coordinate frame. 
    p = thisMagBox.R_SA * ( p - thisMagBox.v_AS );
    % The rest of this function then proceeds with a focus on the source
    % coordinate system until the very end, when the result is rotated back
    % into the "analysis" frame orientation.
end

% Make sure the input position grid looks right.
if size( p, 1 ) ~= 3
    error( 'Input position vector array requires 3 rows' );
end

% Compute integration limit vectors.
rx1 = p(1,:) - min(thisMagBox.vx);
rx2 = p(1,:) - max(thisMagBox.vx);
ry1 = p(2,:) - min(thisMagBox.vy);
ry2 = p(2,:) - max(thisMagBox.vy);
rz1 = p(3,:) - min(thisMagBox.vz);
rz2 = p(3,:) - max(thisMagBox.vz);

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
insidebox = p(1,:)>=min(thisMagBox.vx) & p(1,:)<=max(thisMagBox.vx) & ...
    p(2,:)>=min(thisMagBox.vy) & p(2,:)<=max(thisMagBox.vy) & ...
    p(3,:)>=min(thisMagBox.vz) & p(3,:)<=max(thisMagBox.vz);
Q(:,:,insidebox) = NaN;

% If the user wants the output in the "analysis" frame, we have to rotate
% the resulting Q.
if strcmp(args.Qframe,'A')
    Q = pagemtimes( thisMagBox.R_AS, pagemtimes( Q, thisMagBox.R_SA ) );
end