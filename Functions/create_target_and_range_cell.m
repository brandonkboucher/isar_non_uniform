function [sc,target_locations, ...
    grid, ...
    theta, ...
    u0,...
    is_doppler_aliasing, ...
    Wx] ...
        = create_target_and_range_cell(...
        sc, ...         % scenario parameters
        t_m, ...        % [s] (M x 1) slow-time
        fc, ...         % [Hz] center frequency
        prf, ...        % [Hz] pulse repetition frequency
        f_hat_l, ...    % [Hz] (L x 1) range-frequency
        range ...       % [m] select the range cell
        )

    const = Constants();
    c = const.c;

    %% target definition

    % Define the target's motion
    if sc.is_target_maneuvering ...
            && sc.is_target_accelerating

        w0 = pi; w1 = 1e3; w2 = 1e3; jerk_mag = 5*pi;
        [theta, w0, w1, w2] = create_complex_target_trajectory(...
            0, w0, w1, w2, t_m);

    else
        w0 = pi; % [rad/s] yawing rate
        if sc.is_target_accelerating
            w1 = 1e3; % [rad/s/s] yawing acceleration
            w2 = 1e3;    % [rad/s/s] yawing jerk
        else
            w1 = 0; % [rad/s/s] yawing acceleration
            w2 = 0;    % [rad/s/s] yawing jerk
        end
    
        
        % determine the targets rotation as a function of slow time
        theta = w0 * t_m ...
            + (1/2) * w1 * t_m .* t_m ...
            + (1/3) * w2 * t_m .* t_m .* t_m;
    end

    % determine the angular span which determines the cross
    % range resolution
    sin_span = max(sin(theta)) - min(sin(theta));
    cross_range_resolution = const.c ...
        / (2 * (fc + max(f_hat_l)) * sin_span);
    
    % define the crossrange unambiguous extent, this is the extent at which the 
    % crossrange is unambiguous. if a scatterer's crossrange exceeds this extent
    % then it is ambiguous
    Wx = c * prf / (2 * fc * mean(w0));
    grid.Wx = Wx;

    % determine the number of pixels within the crossrange unambiguous extent
    n_pix_per_amb = round(Wx / cross_range_resolution);
    grid.n_pix_per_amb = n_pix_per_amb;

    % extract the number of ambiguous scatterers 
    n_amb_scat = sc.num_amb_having_scatterers;
    n_amb_if   = sc.num_amb_in_image_former;

    if sc.is_grid_oversampled
        Kc = sc.oversampling_factor;
        grid.sampling_factor = sc.oversampling_factor;
    else
        Kc = 1;
        grid.sampling_factor = 1;
    end

    % redefine the grid in order to handle the number of ambiguous scatterers
    Ny = 1;
    Nx = Kc * n_amb_if * n_pix_per_amb;

    % ensure the number of pixels is odd in both directions
    Nx = Nx + (mod(Nx,2)==0);
    grid.Nx = Nx; grid.Ny = Ny;

    % define the resolution of the grid
    cross_range_pixel_res   = cross_range_resolution / Kc;
    grid.range_pixel_res = 1;
    grid.cross_range_pixel_res = cross_range_pixel_res;
    
    % distance from the origin to the target center
    u0 = 1000; % [m]
    
    % define the range and crossrange grid of the latent image
    % centered on 0 so the critical grid is a subset of the oversampled grid
    y_array = range;
    x_array = ((0:Nx-1) - (Nx-1)/2) * cross_range_pixel_res;
    grid.x_array = x_array; grid.y_array = y_array;
    [X,Y] = meshgrid(x_array,y_array);
    
    % define the x and y component for each grid point
    grid.yk = reshape(Y, [], 1);
    grid.xk = reshape(X, [], 1);
   
    %% sensing matrix and latent image
    
    % number of point scatterers
    Ks = sc.num_of_scatterers;
    
    % assign scatterers depending on the number of ambiguities. if the number
    % of ambiguities is 1 -> [0], if 3 -> [-1,0,1], etc
    ii        = 0:(n_amb_scat-1);
    amb_index = ceil(ii/2) .* (-1).^ii;
    amb_index = sort(amb_index);

    if sc.is_closely_spaced

        [target_locations, amb_of_k] = create_closely_spaced_target_scatterers(...
            Ks, ...
            amb_index, ...
            Wx, ...
            max(y_array) - min(y_array), ...
            cross_range_pixel_res);

    else

        % round-robin so every requested ambiguity holds at least one scatterer
        amb_of_k = amb_index(mod(0:Ks-1, n_amb_scat) + 1).';
    
        % define the target scatterer locations, in meters. Columns are
        % [crossrange, range] to match compute_atom's (x, y) argument order (as used
        % for xk, yk below). Crossrange is uniform within the scatterer's own
        % ambiguity; range is uniform over the grid extent (range has no ambiguity
        % structure here -- it is set by bandwidth, not by the PRF).

        target_locations = [ ...
            rand(Ks,1) * Wx - Wx/2 + amb_of_k * Wx, ...
            (rand(Ks,1) - 0.5) * (y_array(end) - y_array(1)) ];
    end

    % On-Grid rows additionally snap onto the critically-sampled grid. The
    % ambiguity assignment above already happened, so on-grid scatterers are
    % distributed across ambiguities exactly like the off-grid ones.
    if ~sc.is_off_grid
        target_locations(:,1) = round(target_locations(:,1)/cross_range_resolution) * cross_range_resolution;
        target_locations(:,2) = range;
    end

    % The latent image is supported only on the grid, which spans the first
    % n_amb_if ambiguities (same ordering as the scatterer assignment above).
    % Scatterers outside those ambiguities have no atom that represents them:
    % their true atoms correlate with the dictionary only at its sidelobe floor
    % and are not sparsely representable, so they act as unmodeled interference
    % rather than as ghosts folded into ambiguity 0. Reconstruction is scored
    % against the representable scatterers alone.
    jj               = 0:(n_amb_if-1);
    amb_if           = ceil(jj/2) .* (-1).^jj;
    is_latent        = ismember(amb_of_k, amb_if);
    % Ks_latent        = sum(is_latent);
    sc.num_of_latent_scatterers        = Ks;
    latent_locations = target_locations(is_latent, :);

    % determine if Doppler aliasing will occur
    if sc.is_target_maneuvering && sc.is_target_accelerating
        theta_dot = w0 + w1.*t_m + w2.*t_m.^2;
    else
        theta_dot = w0;
    end

    fd_max = max((2*fc/c) * abs(target_locations(:,1)) * max(abs(theta_dot)));
    fprintf('W_x = %.3f m (%d critical pixels); scatterers in ambiguities %s (requested %d)\n', ...
        Wx, n_pix_per_amb, mat2str(unique(amb_of_k).'), n_amb_scat);
    fprintf('image former spans %.2f ambiguities (requested %d), covering %s\n', ...
        Nx*cross_range_pixel_res/Wx, n_amb_if, mat2str(sort(amb_if)));
    fprintf('%d of %d scatterers are representable in the latent image; %d act as interference\n', ...
        sc.num_of_latent_scatterers, Ks, Ks - sc.num_of_latent_scatterers);
    is_doppler_aliasing = false;
    if fd_max > prf/2
        disp('Doppler aliasing will occur')
        is_doppler_aliasing = true;
    else
        disp('No Doppler aliasing')
    end

end

