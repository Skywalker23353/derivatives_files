function plot_sensitivity_comparison(comb_data, noz_data, field_name, label, output_dir, figidx)
% Plot comparison of original and interpolated sensitivity fields using myutils.plot_field
% 
% Inputs:
%   comb_data - combustor sensitivity data structure
%   noz_data - nozzle sensitivity data structure  
%   field_name - name of the sensitivity field to plot
%   label - LaTeX label for the sensitivity field
%   output_dir - directory to save plots
%   figidx - figure index (optional, default based on field)

    if ~isfield(comb_data, field_name)
        fprintf('Warning: Field %s not found in combustor data\n', field_name);
        return;
    end
    
    D = 2e-3; % Diameter scale
   
        % Create combined plot using myutils.plot_field format
        figure(figidx);
        
        % Plot combustor sensitivity field
        myutils.plot_field(figidx, comb_data.R1/D, comb_data.Z1/D, comb_data.(field_name), label);
        hold on;

        myutils.plot_field(figidx, noz_data.R2/D, noz_data.Z2/D, noz_data.(field_name), label);
        
        % Apply formatting
%         xlim([0 3]);
        pbaspect([1 2 1]);
        
        
        % Save plot
        if nargin > 4 && ~isempty(output_dir)
            if ~exist(output_dir, 'dir')
                mkdir(output_dir);
            end
%             saveas(gcf, fullfile(output_dir, sprintf('%s_sensitivity_combined.png', field_name)));
%             saveas(gcf, fullfile(output_dir, sprintf('%s_sensitivity_combined.fig', field_name)));
        end
        
end
