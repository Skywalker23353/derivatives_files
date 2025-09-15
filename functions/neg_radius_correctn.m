function [mat] = neg_radius_correctn(r_mat)
    [rc,cc,pc] =size(r_mat);
    mat = r_mat;
    for r=1:rc
        for c=1:cc
            for p=1:pc
                if (r_mat(r,c,p) < 0)
                    mat(r,c,p) = -1*r_mat(r,c,p);
                else
                    mat(r,c,p) = r_mat(r,c,p);
                end
            end
        end
    end
    clear r_mat;
end