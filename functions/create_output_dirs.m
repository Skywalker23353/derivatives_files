function dirs = create_output_dirs(work_dir)
    % Create all output directories and return their paths
    dirs.sensitivities_dir = fullfile(work_dir, 'species_sensitivities');
    dirs.figures_dir = fullfile(work_dir, 'species_figures');
    dirs.surface_figures_dir = fullfile(work_dir, 'surface_figures');

    % Create directories if they don't exist
    dir_names = fieldnames(dirs);
    for i = 1:length(dir_names)
        dir_path = dirs.(dir_names{i});
        if ~exist(dir_path, 'dir')
            mkdir(dir_path);
            fprintf('Created %s directory: %s\n', dir_names{i}, dir_path);
        end
    end
end