function [fields_struct, successful_fields] = compute_field_derivatives(field_configs, data_dir, C_MAT, Z_MAT, variable_ref_val_list )
    % Compute derivatives for all field data (denominators)
    fprintf('Computing field derivatives (denominators)...\n');
    fields_struct = struct();
    successful_fields = {};
    
    rho_ref = variable_ref_val_list{1,2};
    T_ref = variable_ref_val_list{2,2};
    
    for i = 1:size(field_configs, 1)
        field_name = field_configs{i, 1};
        smooth_suffix = field_configs{i, 2};
        latex_label = field_configs{i, 3};
        plot_title = field_configs{i, 4};
        short_name = field_configs{i, 5};

        fprintf('  Loading %s...\n', field_name);

        % Construct field filename
        field_filename = sprintf('%s%s', field_name, smooth_suffix);
        field_path = sprintf('%s/%s.mat', data_dir, field_filename);

        % Check if field exists
        if ~exist(field_path, 'file')
            fprintf('    Field file not found: %s\n', field_path);
            % Try without smooth suffix
            field_filename_alt = field_name;
            field_path_alt = sprintf('%s/%s.mat', data_dir, field_filename_alt);
            if exist(field_path_alt, 'file')
                fprintf('    Found alternative: %s\n', field_filename_alt);
                field_path = field_path_alt;
            else
                fprintf('    Skipping field: %s\n', field_name);
                continue;
            end
        end

        try
            % Load field data
            loaded_data = load(field_path);
            field = compute_dfdr(loaded_data.DF, C_MAT);
            if strcmp(field_name,'Temperature')
                field = field./T_ref;
            elseif strcmp(field_name,'density')
                field = field./rho_ref;
            end
            
            % Store field information in fields_struct
            fields_struct.(field_name).actual_data = loaded_data.DF;
            fields_struct.(field_name).derivative_wrt_C = field;
            fields_struct.(field_name).latex_label = latex_label;
            fields_struct.(field_name).plot_title = plot_title;
            fields_struct.(field_name).short_name = short_name;
            fields_struct.(field_name).C_MAT = C_MAT;
            fields_struct.(field_name).Z_MAT = Z_MAT;

            successful_fields{end+1} = field_name;
            fprintf('    ✓ Successfully loaded %s\n', field_name);

        catch ME
            fprintf('    ✗ Error loading %s: %s\n', field_name, ME.message);
            continue;
        end
    end
end