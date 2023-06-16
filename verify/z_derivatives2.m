function dudz = z_derivatives2(U, x, y, z)

    x_len = length(x);
    y_len = length(y);
    z_len = length(z);
    
    dz = z(2) - z(1);
    dudz = zeros(x_len, y_len, z_len);
    for i_x = 1:x_len
        for i_y = 1:y_len
            for i_z = 2:z_len - 1
                dudz(i_x, i_y, i_z) = (U(i_x,i_y,i_z-1)-2*U(i_x,i_y,i_z)+U(i_x,i_y,i_z+1))/(dz*dz);
                dudz(i_x, i_y, 1) = (U(i_x,i_y,z_len)-2*U(i_x,i_y,1)+U(i_x,i_y,2))/(dz*dz);
                dudz(i_x, i_y, z_len) = (U(i_x, i_y, z_len-1)-2*U(i_x,i_y,z_len)+U(i_x,i_y,1))/(dz*dz);
            end
        end
    end
end