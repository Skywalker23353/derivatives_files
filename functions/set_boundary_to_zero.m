function field_out = set_boundary_to_zero(field_in, varargin)
% SET_BOUNDARY_TO_ZERO - Sets boundary points of a 2D field to zero
%
% SYNTAX:
%   field_out = set_boundary_to_zero(field_in)
%   field_out = set_boundary_to_zero(field_in, 'BoundaryWidth', width)
%   field_out = set_boundary_to_zero(field_in, 'Boundaries', boundaries)
%   field_out = set_boundary_to_zero(field_in, 'BoundaryWidth', width, 'Boundaries', boundaries)
%
% INPUTS:
%   field_in      - 2D matrix representing the field data
%   
% OPTIONAL PARAMETERS (Name-Value pairs):
%   'BoundaryWidth' - Width of boundary region to set to zero (default: 1)
%                     Can be scalar (same for all boundaries) or 
%                     4-element vector [top, bottom, left, right]
%   
%   'Boundaries'    - Cell array specifying which boundaries to modify
%                     Options: 'top', 'bottom', 'left', 'right', 'all'
%                     Default: {'all'}
%                     Examples: {'top', 'bottom'} or {'left'} or {'all'}
%
% OUTPUTS:
%   field_out     - 2D matrix with specified boundary points set to zero
%
% EXAMPLES:
%   % Set all boundaries (1 point wide) to zero
%   field_zero = set_boundary_to_zero(field);
%
%   % Set top and bottom boundaries (2 points wide) to zero
%   field_zero = set_boundary_to_zero(field, 'BoundaryWidth', 2, 'Boundaries', {'top', 'bottom'});
%
%   % Set different widths for different boundaries
%   field_zero = set_boundary_to_zero(field, 'BoundaryWidth', [3, 2, 1, 1]);
%
%   % Set only left boundary to zero
%   field_zero = set_boundary_to_zero(field, 'Boundaries', {'left'});
%
% NOTES:
%   - Function preserves the original field dimensions
%   - If BoundaryWidth is larger than field dimension, entire rows/columns are set to zero
%   - Useful for CFD post-processing to eliminate boundary artifacts
%   - Compatible with MATLAB conditional statistics and LES field processing
%
% See also: smoothen_fields, thresholding_and_smoothing

% Author: Satyam 
% Date: September 2025
% Part of MATLAB CFD Post-processing Suite

%% Input validation and default values
if ~ismatrix(field_in) || ndims(field_in) ~= 2
    error('Input field must be a 2D matrix');
end

if ~isnumeric(field_in)
    error('Input field must be numeric');
end

[rows, cols] = size(field_in);

% Parse optional inputs
p = inputParser;
addParameter(p, 'BoundaryWidth', 1, @(x) isnumeric(x) && all(x > 0) && (isscalar(x) || length(x) == 4));
addParameter(p, 'Boundaries', {'all'}, @(x) iscell(x) || ischar(x) || isstring(x));
parse(p, varargin{:});

boundary_width = p.Results.BoundaryWidth;
boundaries = p.Results.Boundaries;

% Convert boundaries to cell array if needed
if ischar(boundaries) || isstring(boundaries)
    boundaries = {char(boundaries)};
end

% Validate boundary specifications
valid_boundaries = {'top', 'bottom', 'left', 'right', 'all'};
for i = 1:length(boundaries)
    if ~any(strcmpi(boundaries{i}, valid_boundaries))
        error('Invalid boundary specification: %s. Valid options: %s', ...
              boundaries{i}, strjoin(valid_boundaries, ', '));
    end
end

% Handle 'all' boundary specification
if any(strcmpi(boundaries, 'all'))
    boundaries = {'top', 'bottom', 'left', 'right'};
end

%% Process boundary width parameter
if isscalar(boundary_width)
    % Same width for all boundaries
    width_top = boundary_width;
    width_bottom = boundary_width;
    width_left = boundary_width;
    width_right = boundary_width;
else
    % Different widths: [top, bottom, left, right]
    width_top = boundary_width(1);
    width_bottom = boundary_width(2);
    width_left = boundary_width(3);
    width_right = boundary_width(4);
end

% Ensure widths don't exceed field dimensions
width_top = min(width_top, rows);
width_bottom = min(width_bottom, rows);
width_left = min(width_left, cols);
width_right = min(width_right, cols);

%% Apply boundary modifications
field_out = field_in;

for i = 1:length(boundaries)
    boundary = lower(boundaries{i});
    
    switch boundary
        case 'top'
            if width_top > 0
                field_out(1:width_top, :) = 0;
                fprintf('Set top boundary (%d rows) to zero\n', width_top);
            end
            
        case 'bottom'
            if width_bottom > 0
                field_out(end-width_bottom+1:end, :) = 0;
                fprintf('Set bottom boundary (%d rows) to zero\n', width_bottom);
            end
            
        case 'left'
            if width_left > 0
                field_out(:, 1:width_left) = 0;
                fprintf('Set left boundary (%d columns) to zero\n', width_left);
            end
            
        case 'right'
            if width_right > 0
                field_out(:, end-width_right+1:end) = 0;
                fprintf('Set right boundary (%d columns) to zero\n', width_right);
            end
            
        otherwise
            warning('Unknown boundary specification: %s', boundary);
    end
end

%% Summary output
if nargout == 0
    fprintf('\nBoundary modification summary:\n');
    fprintf('Original field size: %d x %d\n', rows, cols);
    fprintf('Modified boundaries: %s\n', strjoin(boundaries, ', '));
    
    % Count modified points
    modified_points = sum(sum(field_in ~= 0 & field_out == 0));
    total_points = rows * cols;
    fprintf('Points set to zero: %d (%.1f%%)\n', modified_points, 100*modified_points/total_points);
end

end