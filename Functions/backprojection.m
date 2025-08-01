function [rx_signal_bp, output_struct] = ...
    backprojection( ...
        signal, ...
        rx_signal, ...
        ranges, ...
        range_array, ...
        target, ...
        roll, ...
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
    R_roll = [cos(roll), zeros(size(roll,1),1), sin(roll);
     zeros(size(roll,1),1),ones(size(roll,1),1), zeros(size(roll,1),1);
    -sin(roll), zeros(size(roll,1),1), cos(roll)];

    % reshape to dimensions [nPulses x 3 x 3]
    R_roll = reshape(R_roll, [size(roll,1),3,3]);

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

    % formulate the mesh grid
    [X,Y] = meshgrid(x_array, y_array);
    Z = 0;

    % initialize final image
    rx_signal_bp = zeros(Nx, Ny);

    t = tic;
    % iterate through each pixel
    for ix = 1:Nx
        
        if mod(ix, round(Nx/100)) == 0
            
            fprintf('       %3.1f percent complete: %4.2f seconds\n', (ix*100/round(Nx)), toc(t))
            t = tic;
        end

        for iy = 1:Ny

            % extract the pixel location relative to the
            % center of the target
            pixel_location = [X(iy,ix), Y(iy,ix), Z];

            % initialize the image value as zero
            image_value = 0;

            % iterate through each pulse
            for ipulse = 1:signal.num_pulses

                % calculate the location of the pixel
                % relative to the radar
                R = squeeze(R_roll(ipulse,:,:));
                pixel_location_radar = ...
                    pixel_location * R ...
                    + target_positions(ipulse,:);

                % find the range of the pixel to the radar
                range = norm(pixel_location_radar);

                % interpolate the range
                echo = interp1(range_array, ...
                    rx_signal(ipulse, :), ...
                    range, 'linear', 0);

                % Apply phase correction (match propagation delay)
                phase = exp(-1j * 4 * pi * signal.fc * range / const.c);
    
                % Sum contribution
                image_value = image_value + echo * phase;

            end
            rx_signal_bp(ix,iy) = image_value;

        end
    end

    f = figure;
    imagesc(x_array, y_array, abs(rx_signal_bp))
    title('Backprojection image', 'FontSize', 24)
    xlabel('x (cross-range)', 'FontSize', 16)
    ylabel('y (range)', 'FontSize', 16)
    axis square
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(f, 'plots/backprojection.png')

end

