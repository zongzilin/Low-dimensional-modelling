function dudx = z_derivatives(U, x, y, z)

    x_len = length(x);
    y_len = length(y);
    z_len = length(z);
    
    dz = z(2) - z(1);
    dudx = zeros(x_len, y_len, z_len);
    for i_x = 1:x_len
        for i_y = 1:y_len
            for i_z = 2:z_len - 1
                dudx(i_x, i_y, i_z) = (-U(i_x, i_y, i_z - 1) + U(i_x, i_y, i_z + 1))/(2*dz);
                dudx(i_x, i_y, 1) = (-U(i_x, i_y, z_len) + U(i_x, i_y, 2))/(2*dz);
                dudx(i_x, i_y, z_len) = (-U(i_x, i_y, z_len-1) + U(i_x, i_y, 1))/(2*dz);                
            end
        end
    end
end