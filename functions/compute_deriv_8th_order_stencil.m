function [d] = compute_deriv_8th_order_stencil(u)

    N = length(u);
    d = zeros(1,N);
    delta_eta = 1/(N-1);

    coeff = [-25.0/12.0,  -1.0/4.0,   1.0/12.0,   -1.0/60.0,    1.0/280.0,           0,          0,          0,         0;
                    4.0,  -5.0/6.0,   -2.0/3.0,    3.0/20.0,   -4.0/105.0,           0,          0,          0,         0;
                   -3.0,   3.0/2.0,        0.0,    -3.0/4.0,      1.0/5.0,   -1.0/60.0,          0,          0,         0;
                4.0/3.0,  -1.0/2.0,    2.0/3.0,           0,     -4.0/5.0,    3.0/20.0,          0,          0,         0;
               -1.0/4.0,  1.0/12.0,  -1.0/12.0,     3.0/4.0,          0.0,    -3.0/4.0,   1.0/12.0,  -1.0/12.0,   1.0/4.0;
                      0,         0,          0,   -3.0/20.0,      4.0/5.0,         0.0,   -2.0/3.0,    1.0/2.0,  -4.0/3.0;
                      0,         0,          0,    1.0/60.0,     -1.0/5.0,     3.0/4.0,        0.0,   -3.0/2.0,       3.0;
                      0,         0,          0,           0,    4.0/105.0,   -3.0/20.0,    2.0/3.0,    5.0/6.0,      -4.0;
                      0,         0,          0,           0,   -1.0/280.0,    1.0/60.0,  -1.0/12.0,    1.0/4.0, 25.0/12.0;];
   
    for ii = 1:N   
        if (ii <= (4) || ii >= (N - 3))    
            if (ii == 1)
                field = [u(ii:ii+4)];
                d(ii) = (field*coeff(1:5,1))/delta_eta;
                clear field;
            elseif (ii == 2)
                field = [u(ii-1:ii+3)];
                d(ii) = (field*coeff(1:5,2))/delta_eta;
                clear field;
            elseif (ii == 3)
                field = [u(ii-2:ii+2)];
                d(ii) = (field*coeff(1:5,3))/delta_eta;
                clear field;
            elseif (ii == 4)
                field = [u(ii-3:ii+3)];
                d(ii) = (field*coeff(1:7,4))/delta_eta;
                clear field;
            elseif (ii == N - 3)
                field = [u(ii-3:ii+3)];
                d(ii) = (field*coeff(3:9,6))/delta_eta;
                clear field;
            elseif (ii == N - 2)
                field = [u(ii-2:ii+2)];
                d(ii) = (field*coeff(5:9,7))/delta_eta;
                clear field;
            elseif (ii == N - 1)
                field = [u(ii-3:ii+1)];
                d(ii) = (field*coeff(5:9,8))/delta_eta;
                clear field;
            elseif (ii == N)
                field = [u(ii-4:ii)];
                d(ii) = (field*coeff(5:9,9))/delta_eta;
                clear field;
            else 
                fprintf("Error when checking for x location near walls");
            end
        else%if (ii > 4 || ii < (N - 3))
            field = [u(ii-4:ii+4)];
            d(ii) = (field*coeff(1:9,5))/delta_eta;
            clear field;
%         else
%             fprintf("Unable to locate point anywhere for derivative computation");
        end
    end
end
