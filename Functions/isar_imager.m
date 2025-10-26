

function output = isar_imager(scenario)
    
    %% Preprocssing

    const = Constants();

    simStart = tic;
    
    [rx_signal, ranges, phases,scatterer_positions] = ...
        form_raw_isar_data(...
        scenario.signal,...
        scenario.target);

   % define the cross range resolution for the grid
    %cross_range_resolution =  c / (2 * yaws(end) * fc);
    cross_range_resolution = 1;
    
    % define the image dimensions using the range and
    % cross-range resolution and the extent of the
    % dimensions, 
    cross_range_extent = 20; % [m]
    range_extent = 20; % [m]

    [rx_signal_bp, x_array, y_array] = backprojection(...
        scenario.signal, ...
        scenario.target, ...
        rx_signal, ...
        cross_range_resolution, ...
        cross_range_extent, ...
        range_extent);

    % save to output struct
    output.rx_signal_bp = rx_signal_bp;
    output.x_array = x_array;
    output.y_array = y_array;
    output.ranges = ranges;
    output.scatterer_positions = scatterer_positions;

    %% Plot results
    fprintf('Simulation complete: %4.2f seconds\n', toc(simStart))
    fprintf('Plotting results.\n')
    
end

