function [ Qii, alphas, thetas ] = computeQii( ri, rj, rk )

% Qii = computeQii( ri, rj, rk )
% 
% Compute the value of the Qii function across all N positions described by
% the input arrays. The input arrays must have size 2xN with 
% ri(1,:) > ri(2,:) and so on. The output is a 1x1xN array describing the
% value of Qii at each of the N evaluation points. The reason for the shape
% is to make it easy to combine all the terms of Q into a 3x3xN matrix.
% 
% [ Qii, alphas, thetas ] = computeQii( ri, rj, rk ) provides intermediate
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
thetas = nan(2,2,2,N);
alphas = nan(2,N);

% Perform summation over all 8 vertices. If any of the faces are at
% infinity, then the corresponding angle has to be calculated
% differently to avoid numerical problems. Specifically, if the prism
% is semi-infinite in either the j or k directions, the angle becomes:
% arctan( sign(rjb) * rkc/ria ) where the sign is important because the
% limit of the usual expression as rjb goes to ±infinnity retains the
% sign of rjb. If the prism is semi-infinite in both j and k, the corner
% contributes a ±pi/2, with the sign depending on the signs of rjb and rkc.
% Also, note that we don't have to do anything special if
% the prism is semi-infinite in the i direction because the ria in the
% denominator will just push the corresponding angle towards zero.
% Finally, if ria is zero, the corner contributes a ±pi/2 with the sign
% depending on the signs of rjb and rkc.
for a = [ 1 2 ]
    for b = [ 1 2 ]
        for c = [ 1 2 ]

            % Shorthand for logical indexing.
            ria = ri(a,:);
            rjb = rj(b,:);
            rkc = rk(c,:);
            infia = isinf(ria);
            infjb = isinf(rjb);
            infkc = isinf(rkc);

            % Compute sign corresponding to this corner's angle.
            s_abc = (-1)^(a+b+c);

            % Case 0: We are viewing from the plane of the ia face.
            % Here, each corner contributes a constant of ±pi/2, depending
            % on the signs of rjb and rkc terms.
            riaiszero = (ria==0);
            if any(riaiszero)
                thetas(a,b,c,riaiszero) = s_abc * sign(rjb(riaiszero)) .* sign(rkc(riaiszero)) * pi/2;
            end

            % Case 1: Entirely finite.
            inds1 = ~riaiszero & ~infia & ~infjb & ~infkc;
            if any(inds1)
                Rabc = sqrt( ria.^2 + rjb.^2 + rkc.^2 );
                thetas(a,b,c,inds1) = s_abc * atan( (rjb(inds1).*rkc(inds1)) ./ (ria(inds1).*Rabc(inds1)) );
            end

            % Case 2: The ia face is infinitely far away.
            inds2 = ~riaiszero & infia;
            if any(inds2)
                thetas(a,b,c,inds2) = 0;
            end

            % Case 3.1: Only the jb face is infinitely far away.
            inds31 = ~riaiszero & ~infia & infjb & ~infkc;
            if any(inds31)
                thetas(a,b,c,inds31) = s_abc * sign(rjb(inds31)) .* atan( rkc(inds31)./ria(inds31) );
            end

            % Case 3.2: Only the kc face is infinitely far away.
            inds32 = ~riaiszero & ~infia & ~infjb & infkc;
            if any(inds32)
                thetas(a,b,c,inds32) = s_abc * sign(rkc(inds32)) .* atan( rjb(inds32)./ria(inds32) );
            end

            % Case 3.3: The jb and kc faces are both infinitely far away.
            % This is the quarter-plane limit, which is ±pi/2, depending on
            % the signs of ria, rjb, rkc terms.
            inds33 = ~riaiszero & ~infia & infjb & infkc;
            if any(inds33)
                thetas(a,b,c,inds33) = s_abc * sign(rjb(inds33).*rkc(inds33)./ria(inds33)) * pi/2;
            end

        end
    end
end

% Sum across thetas to get alphas.
alphas(1,:) = sum(thetas(1,:,:,:),[ 2 3 ]);
alphas(2,:) = sum(thetas(2,:,:,:),[ 2 3 ]);

% Sum across alphas to get Qii.
Qii = sum( alphas, 1 );

% This will be a 1xN array. But we want 1x1xN (since the
% full Q is 3x3xN). To fix this, we can permute the dimensions.
Qii = permute( Qii, [ 1 3 2 ] );