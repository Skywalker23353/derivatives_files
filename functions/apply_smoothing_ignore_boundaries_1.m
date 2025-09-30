function smoothed_data = apply_smoothing_ignore_boundaries_1(data, window, ignore_boundaries, boundary_width)
% Apply smoothing while ignoring specified boundaries

field = data;
[rows, cols] = size(field);

% Convert ignore_boundaries to cell array if needed
if ischar(ignore_boundaries) || isstring(ignore_boundaries)
    if strcmpi(ignore_boundaries, 'all')
        ignore_boundaries = {'top', 'bottom', 'left', 'right'};
    else
        ignore_boundaries = {char(ignore_boundaries)};
    end
end

% Process boundary width parameter
if isscalar(boundary_width)
    width_top = boundary_width;
    width_bottom = boundary_width;
    width_left = boundary_width;
    width_right = boundary_width;
else
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

% Determine which columns/rows to ignore
ignore_left = any(strcmpi(ignore_boundaries, 'left'));
ignore_right = any(strcmpi(ignore_boundaries, 'right'));
ignore_top = any(strcmpi(ignore_boundaries, 'top'));
ignore_bottom = any(strcmpi(ignore_boundaries, 'bottom'));

% Extract interior region
start_row = ignore_top * width_top + 1;
end_row = rows - ignore_bottom * width_bottom;
start_col = ignore_left * width_left + 1;
end_col = cols - ignore_right * width_right;

if start_row > end_row || start_col > end_col
    warning('Boundary ignoring would remove entire field. Applying standard smoothing.');
    smoothed_data = apply_smoothing(data, window);
    return;
end

interior = field(start_row:end_row, start_col:end_col);

% Create data structure for interior
interior_data = interior;

% Apply smoothing to interior
smoothed_interior_data = apply_smoothing_1(interior_data, window);

% Reconstruct full field
smoothed_field = field;

% Put smoothed interior back
smoothed_field(start_row:end_row, start_col:end_col) = smoothed_interior_data;

smoothed_data = smoothed_field;

fprintf('  Ignored boundaries: %s (width: %s)\n', strjoin(ignore_boundaries, ', '), mat2str(boundary_width));
end