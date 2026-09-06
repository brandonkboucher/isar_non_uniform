function [sc,target_locations, ...
    grid, ...
    theta, ...
    u0,...
    is_doppler_aliasing] ...
        = create_target_and_grid(...
        sc, ...         % scenario parameters
        t_m, ...        % [s] (M x 1) slow-time
        fc, ...         % [Hz] center frequency
        prf, ...        % [Hz] pulse repetition frequency
        f_hat_l, ...    % [Hz] (L x 1) range-frequency
        N_critical ...  % dimension in x and y (for critically sampled)
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
            % Yaw acceleration and jerk are what make an aliased atom
            % distinguishable from its doppelganger: shifting crossrange by Wx
            % leaves the residual phase
            %
            %   2*pi*[ w1*m^2/(2*w0*prf) + w2*m^3/(3*w0*prf^2) ],
            %
            % since the w0 term is exactly the 2*pi*m that creates the alias.
            % A half-cycle residual (enough to break the ambiguity without
            % distorting the aperture) needs
            %
            %   w1 = w0*prf/(M-1)^2   or   w2 = 3*w0*prf^2/(2*(M-1)^3).
            %
            % Both scale with prf, so they belong to the scenario rather than
            % being constants here. Beyond ~4x the half-cycle value the ghost
            % coherence saturates and only the imaging geometry degrades.
            if isfield(sc, 'yaw_acceleration')
                w1 = sc.yaw_acceleration; % [rad/s/s] yawing acceleration
            else
                w1 = 1e3; % [rad/s/s] yawing acceleration
            end
            if isfield(sc, 'yaw_jerk')
                w2 = sc.yaw_jerk; % [rad/s/s/s] yawing jerk
            else
                w2 = 1e3; % [rad/s/s/s] yawing jerk
            end
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
    range_resolution = ...
        const.c / (2 * (max(f_hat_l) - min(f_hat_l))); % [m]

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
        Kr = sc.oversampling_factor; Kc = sc.oversampling_factor;
        grid.sampling_factor = sc.oversampling_factor;
    else
        Kr = 1; Kc = 1;
        grid.sampling_factor = 1;
    end

    % redefine the grid in order to handle the number of ambiguous scatterers
    Ny = Kr * N_critical;
    Nx = Kc * n_amb_if * n_pix_per_amb;

    % ensure the number of pixels is odd in both directions
    Nx = Nx + (mod(Nx,2)==0);
    Ny = Ny + (mod(Ny,2)==0);
    grid.Nx = Nx; grid.Ny = Ny;

    % define the resolution of the grid
    range_pixel_res         = range_resolution / Kr;
    cross_range_pixel_res   = cross_range_resolution / Kc;
    grid.range_pixel_res = range_pixel_res;
    grid.cross_range_pixel_res = cross_range_pixel_res;
    
    % distance from the origin to the target center
    u0 = 1000; % [m]
    
    % define the range and crossrange grid of the latent image
    % centered on 0 so the critical grid is a subset of the oversampled grid
    y_array = ((0:Ny-1) - (Ny-1)/2) * range_pixel_res;
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
            amb_of_k * Wx + (rand(Ks,1) - 0.5) * Wx, ...
            (rand(Ks,1) - 0.5) * (y_array(end) - y_array(1)) ];
    end

    % On-Grid rows additionally snap onto the critically-sampled grid. The
    % ambiguity assignment above already happened, so on-grid scatterers are
    % distributed across ambiguities exactly like the off-grid ones.
    if ~sc.is_off_grid
        target_locations(:,1) = round(target_locations(:,1)/cross_range_resolution) * cross_range_resolution;
        target_locations(:,2) = round(target_locations(:,2)/range_resolution)       * range_resolution;
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
    
    fprintf('Targets located at: \n')
    disp(target_locations)

    is_doppler_aliasing = false;
    if fd_max > prf/2
        disp('Doppler aliasing will occur')
        is_doppler_aliasing = true;
    else
        disp('No Doppler aliasing')
    end

    % --- is the range grid inside its own unambiguous window? ------------
    % The range-frequency spacing df sets an unambiguous range c/(2*df), the
    % range-domain counterpart of Wx. A taller grid folds in range.
    df_bin = f_hat_l(2) - f_hat_l(1);
    Ru = c / (2 * df_bin);
    grid.Ru = Ru;
    y_extent = max(y_array) - min(y_array);
    fprintf('unambiguous range window %.2f m; range grid spans %.2f m (%.2f windows)\n', ...
        Ru, y_extent, y_extent/Ru);
    if y_extent > Ru
        warning(['the range grid is taller than the unambiguous range window, ' ...
            'so scatterers will alias in range as well as crossrange. ' ...
            'Reduce N_critical to at most L = numel(f_hat_l).']);
    end

    % --- is that ambiguity genuine, or smeared away? ----------------------
    % Aliasing requires the atom Wx away in crossrange to reproduce the
    % scatterer's own atom. A wide fractional bandwidth breaks it (each range
    % bin folds at its own Wx(f) = c*prf/(2*(fc+f)*w0)), and so does a target
    % that is not compact in range relative to u0 (rotation-induced range
    % migration). Measure the ghost coherence rather than assume it.
    Wx_per_bin = c * prf ./ (2 * (fc + f_hat_l) * mean(w0));
    grid.Wx_per_bin = Wx_per_bin;

    % Scan for the ghost rather than assuming it sits at exactly x + Wx. Yaw
    % acceleration does not destroy the alias, it moves it: the effective fold
    % distance is set by the mean rotation rate over the aperture, so the ghost
    % drifts from Wx as w1 grows (measured 0.2 pixels at w1 = 0 and 4.3 pixels
    % at w1 = 170 rad/s/s). The peak coherence is what governs whether the
    % ambiguity can be resolved; its displacement is what the ambiguity check's
    % search window (options.amb_refine_pixels) has to cover.
    scan = (-8:0.25:8) * cross_range_pixel_res;
    ghost_coherence = zeros(Ks,1);
    ghost_offset_px = zeros(Ks,1);
    for k = 1:Ks
        ak = compute_atom(target_locations(k,1), target_locations(k,2), ...
            u0, theta, f_hat_l, fc, false);
        ak = ak / norm(ak);
        cg = zeros(numel(scan),1);
        for is = 1:numel(scan)
            gk = compute_atom(target_locations(k,1) + Wx + scan(is), ...
                target_locations(k,2), u0, theta, f_hat_l, fc, false);
            cg(is) = abs(gk' * ak) / norm(gk);
        end
        [ghost_coherence(k), ipk] = max(cg);
        ghost_offset_px(k) = scan(ipk) / cross_range_pixel_res;
    end
    grid.ghost_coherence = ghost_coherence;
    grid.ghost_offset_px = ghost_offset_px;

    fprintf('ghost peak coherence per scatterer: %s\n', ...
        num2str(ghost_coherence.', '%.3f  '));
    fprintf(['ghost displaced from +Wx by: %s pixels ' ...
        '-> set options.amb_refine_pixels >= %d\n'], ...
        num2str(ghost_offset_px.', '%+.2f  '), ...
        max(1, ceil(max(abs(ghost_offset_px)) + 1)));

    fprintf(['fractional bandwidth %.3f -> fold distance %.2f to %.2f m ' ...
        '(%.1f crossrange cells of spread)\n'], ...
        (max(f_hat_l) - min(f_hat_l))/fc, min(Wx_per_bin), max(Wx_per_bin), ...
        (max(Wx_per_bin) - min(Wx_per_bin)) / cross_range_resolution);
    if any(ghost_coherence < 0.8)
        warning(['the Doppler ambiguity is broken for at least one scatterer, ' ...
            'so it will smear rather than alias. Reduce the fractional ' ...
            'bandwidth and/or keep the target compact in range relative to u0.']);
    end

end

