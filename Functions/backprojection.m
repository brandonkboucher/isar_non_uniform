

function [rx_signal_bp, x_array, y_array] = ...
    backprojection( ...
        signal, ...
        target, ...
        rx_signal, ...
        image_parameters, ...
        varargin)

    use_estimation = false;
    if ~isempty(varargin)
        use_estimation = true;
        drift_estimation = varargin{1,1};
    end

    const = Constants();


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

    if strcmp(image_parameters.imager, 'weighted')

        % normalized angular sampling weighting based off of
        % Tang_2016 and Zeng_2012, padded first value
        weights = zeros(size(target.yaws,1),1);
        weights(1) = target.yaws(2) - target.yaws(1);
        weights(end) = target.yaws(end) - target.yaws(end-1);
    
        for k = 2:size(target.yaws,1)-1
            weights(k) = (target.yaws(k+1) - target.yaws(k-1))/2;
        end
        weights = weights ./ sum(weights);

    else
        weights = ones(size(target.yaws,1),1);

    end

    % initialize final image
    rx_signal_bp = zeros(Nx, Ny);
    range_matrix = zeros(Nx, Ny, signal.num_pulses);

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
                range_matrix(ix,iy,ipulse) = range;
    
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
                image_value = image_value + weights(ipulse) * pulse_contribution;
                
            end
            rx_signal_bp(ix,iy) = image_value; % [cross range x range]
        end
    end
    % save('data/range_matrix.mat', "range_matrix");
end

