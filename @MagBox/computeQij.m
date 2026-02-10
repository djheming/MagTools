function [ Qij, f_ab, rho ] = computeQij( ri, rj, rk )

% Qij = computeQij( ri, rj, rk )
% 
% Compute the value of the Qij function across all N positions described by
% the input arrays. The input arrays must have size 2xN with 
% ri(1,:) > ri(2,:) and so on. The output is a 1x1xN array describing the
% value of Qij at each of the N evaluation points. The reason for the shape
% is to make it easy to combine all the terms of Q into a 3x3xN matrix.
% 
% [ Qij, f_ab, f_c, ts ] = computeQij( ri, rj, rk ) provides intermediate
% results that may be helpful for debugging.
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
%   Updated: 2025-05-26
%   Doug Hemingway (douglas.hemingway@utexas.edu)
%   University of Texas at Austin
%

% Initialize.
N = size(ri,2);
rho = nan(2,2,N);
f_ab = nan(2,2,N);
Qij_ab = nan(2,2,N);

% Shorthand for clarity and to aid with logical indexing.
rk1 = rk(1,:);
rk2 = rk(2,:);
infk1 = isinf(rk1);
infk2 = isinf(rk2);
infk = infk1 | infk2;
isAbove = rk1>0 & rk2>0;
isBelow = rk1<0 & rk2<0;

% Compute the lateral (ij) hypotenuese terms for each of the four edges.
s = nan(2,2); % For keeping track of the signs of each edge term.
for a = [ 1 2 ]
    for b = [ 1 2 ]
        rho(a,b,:) = hypot( ri(a,:), rj(b,:) );      
        s(a,b) = (-1)^(a+b);
    end
end

% Loop over the four ij edges (pairing up the terms along the k dimension).
for a = [ 1 2 ]
    for b = [ 1 2 ]

        % Rather than work on each corner in isolation, we pair the
        % corners that are separated in the k-direction. That is, we will
        % go to work on each of the four edges, which are index by (a,b).
        % We do this to enable better handling of singularities where the
        % evaluation point is collinear with one of the edges or where the
        % prism has infinite or semi-infinite extent. 

        % Initialize this edge term.
        edgeTerm = zeros(1,N);

        % Extract the relevant rho and reshape for convenience in logical indexing.
        rho_ab = reshape(rho(a,b,:),1,N);
        rhoiszero = (rho_ab==0); % Points that are collinear with edge ab.
        rhoisinf = isinf(rho_ab); % This edge is infinitely far away.
        okrho = ~rhoiszero & ~rhoisinf; % Finite rho.

        % Branch 0: Evaluation points with infinite rho. This occurs when
        % the prism is semi-infinite in either i or j. In these cases,
        % edge ab contributes nothing to the field and so we do nothing
        % (because the edgeTerm has already been initialized to zero). 
        % One note of caution: In cases where the prism has semi-infinite
        % or infinite extent in the k-direction, the edgeTerms will contain
        % infinities that have been dropped from the Branch 1 calculations
        % below. We can only drop them because they generally cancel due to
        % the alternating signs in the summation over (a,b). The one
        % exception to this is when the prism is semi-infinite in BOTH 
        % i and j, in which case three of the four rho(a,b) terms will be
        % zero, leaving only one edge to contribute to Qij. If this edge is
        % missing a top or bottom, there will be a dangling infinity that
        % does not cancel out. We will handle this special case at the very
        % end (see below).

        % Branch 1: Evaluation points with finite rho.
        if any(okrho)

            % Case 1.0: Edge is entirely finite. 
            inds_1_0 = okrho & ~infk1 & ~infk2;
            edgeTerm(inds_1_0) = asinh(rk1(inds_1_0)./rho_ab(inds_1_0)) - asinh(rk2(inds_1_0)./rho_ab(inds_1_0));

            % Case 1.1: Semi-infinite with rk1=+inf, rk2 finite. The prism has no bottom.
            inds_1_1 = okrho & (infk1 & rk1>0) & ~infk2;
            edgeTerm(inds_1_1) = - asinh(rk2(inds_1_1)./rho_ab(inds_1_1)) - log(rho_ab(inds_1_1));

            % Case 1.2: Semi-infinite with rk1 finite, rk2=-inf. The prism has no top.
            inds_1_2 = okrho & ~infk1 & (infk2 & rk2<0);
            edgeTerm(inds_1_2) = asinh(rk1(inds_1_2)./rho_ab(inds_1_2)) - log(rho_ab(inds_1_2));

            % Case 1.3: Edge is infinite (rk1=+inf, rk2=-inf). Prism has no
            % top or bottom. This is the 2D case.
            inds_1_3 = okrho & (infk1 & rk1>0) & (infk2 & rk2<0);
            edgeTerm(inds_1_3) = -2 * log(rho_ab(inds_1_3));

        end

        % Branch 2: Evaluation points where rho-->0 (singularities).
        if any(rhoiszero)

            % Case 2.0: We are touching the edge.
            inds_2_0 = rhoiszero & ~isAbove & ~isBelow;
            edgeTerm(inds_2_0) = Inf;

            % Case 2.1: Edge is entirely finite and we are above or below
            % (not touching).
            inds_2_1 = rhoiszero & ~infk & ( isAbove | isBelow );
            edgeTerm(inds_2_1) = sign(rk1(inds_2_1)) .* log(rk1(inds_2_1)./rk2(inds_2_1));

            % Case 2.2: Prism has no bottom and we are above.
            inds_2_2 = rhoiszero & infk1 & isAbove;
            edgeTerm(inds_2_2) = - log(2*abs(rk2(inds_2_2)));

            % Case 2.3: Prism has no top and we are below.
            inds_2_3 = rhoiszero & infk2 & isBelow;
            edgeTerm(inds_2_3) = - log(2*abs(rk1(inds_2_3)));

        end

        % Apply the approriate sign for this term.
        Qij_ab(a,b,:) = s(a,b) * edgeTerm;
        f_ab(a,b,:) = exp(edgeTerm);

    end

end

% Sum the terms across a and b.
Qij = sum( Qij_ab, [ 1 2 ] );

% Final check: make sure we are not violating our assumption about the rk1
% and rk2 infinities cancelling. They will generally cancel because we are
% summing (with alternating signs) across the four edges. However, if we
% have a prism that is semi-infinite in both the i and j directions, then
% three of our four edges will be contributing nothing to the summation for
% Qij. This means that if the remaining edge (the only contributing edge)
% contains an infinity, there is nothing to cancel it out. We only need to
% look at points where rk1 or rk2 are infinite.
if any(infk)

    % Here, we have a prism that goes to infinity in rk1 and/or rk2. 
    % Now we just need to see how many of the corresponding rhos are also
    % infinite. If it's an even number (0 or 2), we're good because the
    % corresponding infinite rkc terms will cancel. On the other hand, if
    % there is an odd number of rhos that have gone to infinity, we'll be
    % left with a dangling rkc infinity that doesn't cancel. 
    
    % Track the edges that have finite rhos and note their signs.
    sfrho = s .* ~isinf(rho); 

    % For each evaluation point, the sum of sfrho will be zero if all the
    % infinities cancel. If it's not zero, its sign tells us the sign of the
    % infinity (i.e., if the rho(a,b) term is finite, the infinity is
    % positive for even a+b and negative for odd a+b).
    inf_sign = reshape(sum(sfrho,[1 2]),1,N);
    badinf = infk & (inf_sign~=0);
    if any(badinf)
        Qij(badinf) = inf_sign(badinf) * Inf;
    end
    
end





