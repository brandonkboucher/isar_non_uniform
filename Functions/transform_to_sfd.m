function output = transform_to_sfd(output, image_parameters)
    
    % determine the size of the grid
    Nx = length(output.x_array);
    Ny = length(output.y_array);

    % calculate the kx and ky grid points
    output.kx = (-Nx/2:Nx/2-1)...
        /(Nx*image_parameters.x_pixel_resolution);  % m^-1
    output.ky = (-Ny/2:Ny/2-1)....
        /(Ny*image_parameters.y_pixel_resolution);  % m^-1

    % transform from the image domain to k-space
    output.rx_signal_sfd = sd_to_sfd(output.rx_signal_bp);

end

