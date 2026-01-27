function [ Qij, f_ab, f_c, ts, Rsq ] = computeQij( ri, rj, rk )

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
ts = nan(2,2,2,N);
f_ab = nan(2,2,N);
f_c = nan(2,N);
Rsq = nan(2,2,N);

% Perform summation over all 8 vertices. For a finite prism, we
% can then calculate the f_ab ratios and then Qij from there.
% If any of the i or j faces are at infinity, then we can't
% calculate all the f_ab ratios from the ts and instead have to
% calculate some of them directly. If any of the k faces are at
% infinity, we can't use the f_ab ratios and we have to
% calculate the f_c ratios instead.
if ~any( isinf(rk(:,1)) )
    
    % Here, the prism is fully finite in the k-direction. If
    % any of the i or j faces is at infinity, we already know
    % some of the corresponding f_ab ratios.
    for a = [ 1 2 ]
        if isinf(ri(a,1))
            f_ab(a,:,:) = 1;
        else
            for b = [ 1 2 ]
                if isinf(rj(b,1))
                    f_ab(:,b,:) = 1;
                else
                    for c = [ 1 2 ]
                        Rabc = sqrt( ri(a,:).^2 + rj(b,:).^2 + rk(c,:).^2 );
                        ts(a,b,c,:) = rk(c,:) + Rabc;
                    end
                    f_ab(a,b,:) = ts(a,b,1,:) ./ ts(a,b,2,:);
                end
            end
        end
    end

    % We could compute Qij from the f_ab ratios with
    % Qij = log( (f_ab(1,1,:).*f_ab(2,2,:))./(f_ab(1,2,:).*f_ab(2,1,:)) ); 
    % but this does not handle singularities properly. Instead, we'll start
    % from the ts terms. When the evaluation point is aligned with the ij
    % edge of the prism, and if rkc<0, the corresponding ts is zero and so
    % its logarithm will be -Infinity. Outside the prism, these always come
    % in cancelling pairs, so we can simply replace those log(ts) terms
    % with zero.
    logts = log(ts);
    ts0inds = isinf(logts);
    logts(ts0inds) = 0;
    [as, bs, cs] = ndgrid(1:2, 1:2, 1:2);
    S = repmat( -(-1).^(as + bs + cs), 1, 1, 1, N ); % Places correct sign in front of logarithm terms.
    A15terms = S .* logts; % A15 refers to the equation in the paper.
    Qij = sum( A15terms, [ 1 2 3 ] ); 
    
    % Reshape to get the desired 1x1xN layout.
    Qij = reshape( Qij, 1, 1, [] );

else

    % Here, the prism has one or both of its k-faces at
    % infinity. We'll need to use f_c here instead of the
    % f_ab ratios.
    for c = [ 1 2 ]
        if ~isinf(rk(c,1))
            % This is the finite k-face. Calcluate the ts in
            % the normal way and use this to make f_c.
            for a = [ 1 2 ]
                for b = [ 1 2 ]
                    Rabc = sqrt( ri(a,:).^2 + rj(b,:).^2 + rk(c,:).^2 );
                    ts(a,b,c,:) = rk(c,:) + Rabc;
                end
            end
            f_c(c,:) = ( ts(1,1,c,:).*ts(2,2,c,:) )./( ts(1,2,c,:).*ts(2,1,c,:) );
        elseif rk(c,1) > 0
            % This is the case where rk1 goes to +infinity. In that case,
            % all the ts converge to some large value such that f1 = 1.
            f_c(c,:) = 1;
        elseif rk(c,1) < 0
            % This is the case where rk2 goes to -infinity. In that case,
            % each of the ts shrink towards zero (but do not converge) in
            % such a way that f2 approaches a finite number related to the
            % cross-sectional distances to each of the prism's four finite
            % edges.
            for a = [ 1 2 ]
                for b = [ 1 2 ]
                    Rsq(a,b,:) = ri(a,:).^2 + rj(b,:).^2;
                end
            end
            f_c(c,:) = ( Rsq(1,1,:).*Rsq(2,2,:) )./( Rsq(1,2,:).*Rsq(2,1,:) );
        end
    end

    % Compute Qij from the f_c ratios.
    Qij = log( f_c(1,:)./f_c(2,:) );
    
    % This will be a 1xN array. But we want 1x1xN (since the
    % full Q is 3x3xN). To fix this, we can permute the
    % dimensions.
    Qij = permute( Qij, [ 1 3 2 ] );

end