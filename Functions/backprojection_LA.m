

function [rx_signal_bp, x_array, y_array, asm, gim] = ...
    backprojection_LA( ...
        signal, ...
        target, ...
        rx_signal, ...
        image_parameters, ...
        truncate)

    const = Constants();
    c = const.c;

    %% debugging
    [X,Y] = meshgrid(...
        image_parameters.x_array,...
        image_parameters.y_array);
    [trow, tcol] = find(...
         target.scatter_relative_positions(1,1) == X ...
        & target.scatter_relative_positions(1,2) == Y);

    % define a truncated grid around one of the target
    % scatterers
    if image_parameters.truncate
        
        % extract target location (x,y)
        target_location_x = X(trow, tcol);
        target_location_y = Y(trow, tcol);

        % re-formulate x and y array
        x_array = ...
            (target_location_x - image_parameters.truncation_margin)...
            :image_parameters.x_pixel_resolution...
            :(target_location_x + image_parameters.truncation_margin);
        y_array = (target_location_y - image_parameters.truncation_margin)...
            :image_parameters.y_pixel_resolution...
            :(target_location_y + image_parameters.truncation_margin);

        % re-formulate the meshgrid
        [X,Y] = meshgrid(...
            x_array,...
            y_array);

        % redefine x and y array sizes
        Nx = size(x_array, 2);
        Ny = size(y_array, 2);
        
    else

        % set as default, full grid
        Nx = image_parameters.Nx;
        Ny = image_parameters.Ny;
        x_array = image_parameters.x_array;
        y_array = image_parameters.y_array;

    end

    t = tic;

    % iterate through each pixel
    fprintf('Performing backprojection\n')

    % from Fosca 2021, Equation 7: find the indices
    % corresponding to a minimization of the
    % difference between the distance from each grid point to each range
    
    % iterate through each pulse
    asm = sparse(Nx * Ny, Nx * Ny * signal.num_pulses);
    gim_array = cell(signal.num_pulses,1);
    
    d_max = norm(target.target_center_position);
    d_min = d_max;
    
    for ipulse = 1:signal.num_pulses
    
        if mod(ipulse, round(signal.num_pulses/10)) == 0
            fprintf('   %3.1f percent complete: %4.2f seconds\n',...
                (ipulse*100/round(signal.num_pulses)), toc(t))
            t = tic;
        end

        phase = zeros(Nx * Ny, 1);
        gim_per_pulse = zeros(Nx * Ny, signal.num_range_bins);
    
        for ix = 1:Nx % cross range
        
            for iy = 1:Ny % range
        
                % extract the pixel location relative to the
                % center of the target
                pixel_rel_location = [x_array(ix), y_array(iy), 0];
        
                idx = find(...
                     x_array(ix) == X(:) ...
                    & y_array(iy) == Y(:));
            
                % formulate the rotation matrix about the z axis
                R = [cos(target.yaws(ipulse)), -sin(target.yaws(ipulse)), 0; ...
                     sin(target.yaws(ipulse)), cos(target.yaws(ipulse)), 0; ...
                     0, 0, 1];
            
                pixel_location_radar = ...
                    (R * pixel_rel_location')' ...
                    + target.target_center_position;
            
                % find the range of the pixel to the radar
                d = norm(pixel_location_radar); % eq 5
               
                k1 = find(signal.range_array <= d, 1, 'last'); % eq 7.
                if isempty(k1) || k1 > signal.num_range_bins
                    warning('k1 out of bounds.')
                end
                k2 = k1 + 1;
                
                w1 = 1 - (d - signal.range_array(k1)) ...
                    / (signal.range_array(k2) - signal.range_array(k1));
                w2 = 1 - (signal.range_array(k2) - d) ...
                    / (signal.range_array(k2) - signal.range_array(k1));
                gim_per_pulse(idx, [k1, k2]) = [w1, w2];
            
                phase(idx) = exp(1j * 4 * pi * signal.fc * d / c);
    
                if d_max < d
                    d_max = d;
                end
                if d_min > d
                    d_min = d;
                end
            end
        end
        gim_array{ipulse} = gim_per_pulse;
        asm(:, (ipulse - 1)*Nx*Ny + (1:Nx*Ny)) = ...
            spdiags(phase, 0, Nx*Ny, Nx*Ny);
    end
    
    % truncate the gim to prevent storage constraints
    min_idx = find(d_min >= signal.range_array, 1, "last");
    max_idx = find(d_max <= signal.range_array, 1);
    g = cellfun(@(x) x(:,min_idx:max_idx), gim_array, 'UniformOutput', false);
    
    % truncate ranges
    rx_signal_trunc = rx_signal(:,min_idx:max_idx); % [pulses x ranges]
    rx_signal_trunc_vec = reshape(rx_signal_trunc.', [], 1);
    
    % formulate the interpolation matrix 
    % [Nx * Ny * num_pulses x num_pulses * num_ranges]
    gim = blkdiag(g{:});
    
    % calculate the sparsifying transformation matrix Psi
    psi = asm * gim; % [Nx * Ny x num_pulses * num_ranges]
    
    % apply asm for azimuth compression
    rx_signal_bp = psi * rx_signal_trunc_vec; % eq 10
    rx_signal_bp = reshape(rx_signal_bp, [Nx, Ny]);
    rx_signal_bp = rx_signal_bp.';

end

