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
% sign of rjb. Also, note that we don't have to do anything special if
% the prism is semi-infinite in the i direction because the ria in the
% denominator will just push the corresponding angle towards zero.
for a = [ 1 2 ]
    for b = [ 1 2 ]
        for c = [ 1 2 ]
            if isinf(ri(a,1))
                thetas(a,b,c,:) = 0;
            elseif isinf(rj(b,1))
                thetas(a,b,c,:) = (-1)^(a+b+c) * atan( sign(rj(b,1)) * rk(c,:)./ri(a,:) );
            elseif isinf(rk(c,1))
                thetas(a,b,c,:) = (-1)^(a+b+c) * atan( sign(rk(c,1)) * rj(b,:)./ri(a,:) );
            else
                Rabc = sqrt( ri(a,:).^2 + rj(b,:).^2 + rk(c,:).^2 );
                thetas(a,b,c,:) = (-1)^(a+b+c) * atan2( (rj(b,:).*rk(c,:)), (ri(a,:).*Rabc) );
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