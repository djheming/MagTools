function fh = showBfieldVectors( thisSource, varargin )

% fh = showBfieldVectors( thisSource, ... ) make a figure showing
% the vector B field in a 3D volume.
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
args = BaseTools.argarray2struct( varargin, { 'Q', [], 'view', [], 'z_down', false, 'figsize', 520, 'title', false, 'overlay_vectors', true, 'show_source', true } );
optargs = BaseTools.struct2argarray(args);

% Check leading positional arguments for an axes and a SurveyField. 
for k = 1 : length(args.posArgs)
    if isa( args.posArgs{k}, 'matlab.graphics.axis.Axes' )
        ah = args.posArgs{k};
    elseif isa( args.posArgs{k}, 'SurveyField' )
        survey = args.posArgs{k};
    end
end

% If the survey field hasn't been supplied, make one up.
if ~exist( 'survey', 'var' ) || isempty( survey )
    survey = SurveyField( thisSource ); 
end

% If the user has supplied an array of surveys, handle them recursively.
num_surveys = length(survey);
if num_surveys > 1
    fh = cell(num_surveys,1);
    for k = 1 : num_surveys
        fh{k} = thisSource.showBfieldContours( component, survey(k), varargin{:} );
    end
    return;
end

% Compute B field at all evaluation points.
B = thisSource.computeBfield(survey.p,'Q',args.Q);

% If an Axes handle has not been supplied, make one here.
if ~exist( 'ah', 'var' ) || isempty( ah )
    fh = figure;
    ah = axes('Parent',fh);
    xlabel( ah, 'x' );
    ylabel( ah, 'y' );
    zlabel( ah, 'z' );
    hold( ah, 'on' );
    axis( ah, 'equal' );
    axis( ah, 'tight' );
else
    fh = ancestor(ah,'figure');
end
MagSource.drawBfieldVectors( ah, survey, B, optargs{:} );
