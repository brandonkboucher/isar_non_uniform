

function [rx_signal_bp, x_array, y_array] = ...
    backprojection( ...
        signal, ...
        target, ...
        rx_signal, ...
        cross_range_resolution, ...
        cross_range_extent, ...
        range_extent)

    const = Constants();

    %% backprojection

    % the output of the backprojection algorithm will be 
    % [Nx x Ny] where Nx is cross-range, Ny is range
    Nx = round(cross_range_extent / cross_range_resolution);
    Ny = round(range_extent / signal.range_resolution);
    
    % calculate the x and y values that define the image and
    % the pixel locations relative to the radar
    x_array = (-floor(Nx/2):ceil(Nx/2)-1) * cross_range_resolution;
    y_array = (-floor(Ny/2):ceil(Ny/2)-1) * signal.range_resolution;
    
    % initialize final image
    rx_signal_bp = zeros(Nx, Ny);
    
    t = tic;
    
    % iterate through each pixel
    fprintf('Performing backprojection\n')
    for ix = 1:Nx % cross range
        
        if mod(ix, round(Nx/10)) == 0
            fprintf('   %3.1f percent complete: %4.2f seconds\n',...
                (ix*100/round(Nx)), toc(t))
            t = tic;
        end
    
        for iy = 1:Ny % range
    
            % extract the pixel location relative to the
            % center of the target
            pixel_rel_location = [x_array(ix), y_array(iy), 0];
    
            % initialize the image value as zero
            image_value = 0;
    
            % iterate through each pulse
            for ipulse = 1:signal.num_pulses
    
                % formulate the rotation matrix about the z axis
                R = [cos(target.yaws(ipulse)), -sin(target.yaws(ipulse)), 0; ...
                     sin(target.yaws(ipulse)), cos(target.yaws(ipulse)), 0; ...
                     0, 0, 1];
    
                pixel_location_radar = ...
                    (R * pixel_rel_location')' ...
                    + target.target_center_position;
    
                % find the range of the pixel to the radar
                range = norm(pixel_location_radar);
    
                % interpolate the range
                echo = interp1(...
                    signal.range_array, ...
                    rx_signal(ipulse, :), ...
                    range, ...
                    'linear', 0);
    
                % Apply phase correction (match propagation delay
                phase = exp( 1j * 4 * pi * signal.fc * range / const.c);
    
                % calculation the pulse contribution
                pulse_contribution = echo * phase;
    
                % Sum contribution
                image_value = image_value + pulse_contribution;
    
            end
            rx_signal_bp(ix,iy) = image_value; % [cross range x range]
    
        end
    end
end

