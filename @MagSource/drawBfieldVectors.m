function ah = drawBfieldVectors( ah, survey, B, varargin )

% ah = drawBfieldVectors( ah, survey, B, ... ) make a figure showing
% vectors of the supplied field.  
%
%
%   Disclaimer: This code is provided as-is, has been tested only very
%   informally, and may not always behave as intended. I find it useful for
%   my own work, and I hope you will too, but I make no guarantees as to
%   the accuracy or robustness of this code. This code is also actively
%   under development and future versions may not be backward compatible.
%
%   Created: 2025-05-09
%   Doug Hemingway (douglas.hemingway@utexas.edu)
%   University of Texas at Austin
%


% Read input arguments.
args = BaseTools.argarray2struct( [ { 'clr', [ .6 .6 .8 ], 'view', [], 'figsize', 520, 'overlay', true, 'scale_factor', [] } varargin ] );
if ~args.overlay
    cla( ah );
end

% Scale vectors for display purposes.
if isempty( args.scale_factor )
    args.AutoScale = 'on';
    args.scale_factor = 1;
else
    args.AutoScale = 'off';
end
Bx = B(1,:) * args.scale_factor;
By = B(2,:) * args.scale_factor;
Bz = B(3,:) * args.scale_factor;

% Overlay a grid of vectors.
ah = quiver3( ah, survey.p(1,:), survey.p(2,:), survey.p(3,:), Bx, By, Bz, 'AutoScale', args.AutoScale, 'Color', args.clr, 'LineWidth', 2.0 );


