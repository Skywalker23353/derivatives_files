function smooth_field = thresholding_and_smoothing(field,window,varargin)
    if nargin > 3
        th_f_max = varargin{1};
        th_f_min = varargin{2};
    else
        th_f_max = varargin{1};
        th_f_min = varargin{1};
    end
    if nargin > 4
        smoothing_cycle = varargin{3};
    else
        smoothing_cycle = 2;
    end

    max_val = max(max(field));
    min_val = min(min(field));
    field(field >= th_f_max*max_val) = th_f_max*max_val;
    field(field <= th_f_min*min_val) = th_f_min*min_val;
    smooth_field = field;
    for i = 1:smoothing_cycle
        smooth_field = myutils.f_return_smooth_field(smooth_field,window,'row');
        smooth_field = myutils.f_return_smooth_field(smooth_field,window,'col');
    end
end