function [a,i] = check_NAN(mat)
    a = isnan(mat);
    [a,i] = max(max(max(a)));
end
