function dudx = x_derivatives2(U, x, y, z)

    x_len = length(x);
    y_len = length(y);
    z_len = length(z);

    dx = x(2) - x(1); % x coord distribution is even
    dudx = zeros(x_len, y_len, z_len);
    for i_x = 2:x_len - 1
        for i_y = 1:y_len
            for i_z = 1:z_len
                dudx(i_x, i_y, i_z) = (U(i_x-1, i_y, i_z) - 2*U(i_x, i_y, i_z) + U(i_x+1, i_y, i_z))/(dx*dx);
                dudx(1, i_y, i_z) = (U(x_len, i_y, i_z) - 2*U(1, i_y, i_z) + U(2, i_y, i_z))/(dx*dx);
                dudx(x_len, i_y, i_z) = (U(x_len-1, i_y, i_z) - 2*U(x_len, i_y, i_z) + U(1, i_y, i_z))/(dx*dx);
            end
        end
    end

end