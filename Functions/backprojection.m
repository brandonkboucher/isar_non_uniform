function [rx_signal_bp, output_struct] = ...
    backprojection( ...
        signal, ...
        rx_signal, ...
        ranges, ...
        range_array, ...
        target, ...
        yaw, ...
        target_positions, ...
        output_struct)

    const = Constants();

    % the goal of this script is to perform backprojection
    % where for each pixel, the energy from each pulse is
    % coherently summed. The grid is defined from the
    % reference frame of the center of the target

    % Backprojection functions by coherently summing all of
    % the pulses for each pixel of the image. The image
    % plane is target-fixed meaning that the imaging plane
    % is centered at the target and moves with the target.
    % This means for each pixel, we must consider the
    % translational motion of the target and the rotation of
    % the target in order to define the location of each
    % pixel with respect to the target.

    % define the rotational matrix that serves to rotate the
    % imaging plane with the target

    R_yaw = ...
        [cos(yaw), -sin(yaw), zeros(size(yaw,1),1); ...
         sin(yaw), cos(yaw), zeros(size(yaw,1),1); ...
         zeros(size(yaw,1),1), zeros(size(yaw,1),1), ones(size(yaw,1),1)];

    % R_yaw = [cos(yaw), zeros(size(yaw,1),1), sin(yaw);
    %  zeros(size(yaw,1),1),ones(size(yaw,1),1), zeros(size(yaw,1),1);
    % -sin(yaw), zeros(size(yaw,1),1), cos(yaw)];

    % reshape to dimensions [nPulses x 3 x 3]
    R_yaw = reshape(R_yaw, [size(yaw,1),3,3]);

    % define the image dimensions using the range and
    % cross-range resolution and the extent of the
    % dimensions, 
    cross_range_extent = 100; % [m]
    range_extent = 100; % [m]
    
    % the output of the backprojection algorithm will be 
    % [Nx x Ny] where Nx is cross-range, Ny is range
    Nx = round(cross_range_extent / signal.cross_range_resolution);
    Ny = round(range_extent / signal.range_resolution);

    % calculate the x and y values that define the image and
    % the pixel locations relative to the radar
    x_array = linspace(-cross_range_extent/2, ...
        cross_range_extent/2, Nx)';
    
    y_array = linspace(-range_extent/2, ...
        range_extent/2, Ny)';

    % initialize final image
    rx_signal_bp = zeros(Nx, Ny);

    t = tic;
    
    % f = waitbar(0, 'Performing backprojection');
    % iterate through each pixel
    for ix = 1:Nx % cross range
        
        if mod(ix, round(Nx/100)) == 0
            % waitbar(ix/round(Nx),f,sprintf('Progress: %d %%', floor(ix/Nx*100)))
            fprintf('       %3.1f percent complete: %4.2f seconds\n',...
                (ix*100/round(Nx)), toc(t))
            t = tic;
        end

        for iy = 1:Ny % range

            % extract the pixel location relative to the
            % center of the target
            pixel_location = [x_array(ix), y_array(iy), 0];

            % initialize the image value as zero
            image_value = 0;

            % iterate through each pulse
            for ipulse = 1:signal.num_pulses

                % calculate the location of the pixel
                % relative to the radar
                R = squeeze(R_yaw(ipulse,:,:));
                pixel_location_radar = ...
                    (R * pixel_location')' ...
                    + target_positions(ipulse,:);

                % find the range of the pixel to the radar
                range = norm(pixel_location_radar);

                % interpolate the range
                echo = interp1(...
                    range_array, ...
                    rx_signal(ipulse, :), ...
                    range, ...
                    'linear', ...
                    0);

                % Apply phase correction (match propagation delay)
                phase = ...
                    exp(1j * 4 * pi * signal.fc * range / const.c);
    
                % calculation the pulse contribution
                pulse_contribution = echo * phase;

                % Sum contribution
                image_value = ...
                    image_value + pulse_contribution;

            end
            rx_signal_bp(ix,iy) = image_value; % [cross range x range]

        end
    end

    output_struct.rx_bp = rx_signal_bp;
    output_struct.x_bp = x_array;
    output_struct.y_bp = y_array;

end

