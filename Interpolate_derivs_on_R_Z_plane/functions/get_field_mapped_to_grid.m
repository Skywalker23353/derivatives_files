function field_mat = get_field_mapped_to_grid(NI,NJ,Avg_field)

    field_vec = Avg_field(:,1);
    field_mat = zeros(NI,NJ);
    indx = 1;
    for ii=1:NI
        for jj=1:NJ
            field_mat(ii,jj) = field_vec(indx,1);
            indx = indx + 1;
        end
    end
    clear field_vec;

end
