function field_out = custom_thresholding_and_smoothing(field_in, sensitivities_dir, species_name, field_name, window)

    % Default output (in case of early exit)
    field_out = field_in;

    % Construct parameter file path
    thld_and_smth_parameters_file = fullfile(sensitivities_dir, 'custom_thld_and_smth_param.mat');

    % Check file existence
    if ~exist(thld_and_smth_parameters_file, 'file')
        warning('ERROR: Thresholding and smoothing parameter file not found: %s', thld_and_smth_parameters_file);
        return;
    end

    % Load parameter file
    S = load(thld_and_smth_parameters_file);
    if ~isfield(S, 'data')
        warning(' ERROR: Parameter file does not contain ''data'' variable.');
        return;
    end

    % Build dataset name
    dsetName = sprintf('d%s_d%s', species_name, field_name);

    % Look up parameters
    matchMask = strcmp(S.data(:,1), dsetName);
    if ~any(matchMask)
        warning('No parameters found for dataset: %s', dsetName);
        return;
    end

    % Extract parameter row
    pmt = S.data(matchMask, :);

    % Extract thresholds safely
    maxTh = pmt{2};
    minTh = pmt{3};
    nCycles = pmt{4};

    % Logging (format numbers correctly)
    fprintf('Smoothing parameters for %s:\n', dsetName);
    fprintf('  Max threshold    : %.4f\n', maxTh);
    fprintf('  Min threshold    : %.4f\n', minTh);
    fprintf('  Smoothing cycles : %d\n', nCycles);

    % Apply thresholding and smoothing
    field_out = thresholding_and_smoothing(field_in, window, maxTh, minTh, nCycles);
end
