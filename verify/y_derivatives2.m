function dudy = y_derivatives2(U, x, y, z, dw2)
    x_len = length(x);
    y_len = length(y);
    z_len = length(z);
    
    dudy = zeros(x_len, y_len, z_len);
    for i_x = 1:x_len
        for i_z = 1:z_len
            dudy(i_x, :, i_z) = dw2(:,:,2)*U(i_x,:,i_z)';
        end
    end
end