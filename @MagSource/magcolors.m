function cmap = magcolors( user_N )

% Define a colormap that runs from dark cool colors, to pure white in the
% middle, then to dark warm colors.
%
%   Disclaimer: This code is provided as-is, has been tested only very
%   informally, and may not always behave as intended. I find it useful for
%   my own work, and I hope you will too, but I make no guarantees as to
%   the accuracy or robustness of this code. This code is also actively
%   under development and future versions may not be backward compatible.
%
%   Doug Hemingway (djheming@berkeley.edu)
%   University of California, Berkeley
%   2017-07-31
%

internal_N = 101;
x = linspace( -1, 1, internal_N );
xn = x( x < 0 );
xp = x( x >= 0 );
yn = [ 1-0.9*abs(xn).^0.5; 1-0.8*abs(xn).^1.2; 1-0.6*abs(xn).^2 ];
yp = [ 1-0.4*abs(xp).^3; 1-0.8*abs(xp); 1-0.9*abs(xp).^0.3 ];
cmap = [ yn yp ]';
if exist( 'user_N', 'var' ) && user_N ~= internal_N
    xq = linspace( -1, 1, user_N )';
    cmap = [ interp1(x,cmap(:,1),xq) interp1(x,cmap(:,2),xq) interp1(x,cmap(:,3),xq) ];
end
% figure; plot( cmap );

