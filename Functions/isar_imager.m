

function output = isar_imager( ...
    scenario, ...
    varargin)
    
    simStart = tic;
    const = Constants();

    % calculate and display the range and cross range
    % resolution
    cross_range_resolution = const.c ...
        / (2 * scenario.target.yaws(end) * scenario.signal.fc);
    range_resolution = scenario.signal.range_resolution;

    fprintf('Range resolution:          %2.4f meters\n', range_resolution);
    fprintf('Crossrange resolution:     %2.4f meters\n', cross_range_resolution);
    
    if ~isempty(varargin)
        image_parameters = varargin{1};
    else
        
        image_parameters.x_pixel_resolution = ...
            range_resolution;
        image_parameters.y_pixel_resolution = ...
            cross_range_resolution;

        % define the image dimensions using the range and
        % cross-range resolution and the extent of the
        % dimensions, 
        image_parameters.x_extent = 7; % [m]
        image_parameters.y_extent = 7; % [m]
        
    end

    % define the backprojection grid
    image_parameters = ...
        define_backprojection_grid(image_parameters);

    fprintf('x pixel resolution:                        %2.4f meters\n', image_parameters.x_pixel_resolution);
    fprintf('y pixel resolution:                        %2.4f meters\n', image_parameters.y_pixel_resolution);
    fprintf('Pixels per range resolution cell:          %2.4f\n', image_parameters.x_pixel_resolution / range_resolution);
    fprintf('Pixels per cross-range resolution cell:    %2.4f\n', image_parameters.y_pixel_resolution / cross_range_resolution);
    
    fprintf('-------------------------------------------\n')

    % form the raw image
    [rx_signal, ranges, ~, scatterer_positions] = ...
        form_raw_isar_data(...
        scenario.signal,...
        scenario.target);
    
    [rx_signal_bp, x_array, y_array] = backprojection(...
        scenario.signal, ...
        scenario.target, ...
        rx_signal, ...
        image_parameters, ...
        true);

    % save to output struct
    output.rx_signal = rx_signal;
    output.rx_signal_bp = rx_signal_bp;
    output.x_array = x_array;
    output.y_array = y_array;
    output.ranges = ranges;
    output.scatterer_positions = scatterer_positions;

    output = transform_to_sfd(output, image_parameters);

    %% Plot results
    fprintf('Simulation complete: %4.2f seconds\n', toc(simStart))
    fprintf('Plotting results.\n')
    
end

