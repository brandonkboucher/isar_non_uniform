function output_struct = isar_imager_bp(signal,target)
    
    % instantiate constants
    const = Constants;

    simStart = tic;
    
    %% Define the tx signal model
    
    % calculate the angular extent over the course of the target
    % flyout
    yaw = target.yaw .* signal.t_m ...   
        + (1/2) * target.yawing_rate .* signal.t_m .^ 2;
    angular_extent = yaw(end);
    if angular_extent > 2*pi; angular_extent = 2*pi; end
    signal.cross_range_resolution = ...
        const.c / (2 * angular_extent * signal.fc);
    
    % save data
    output_struct.signal = signal;
    output_struct.target = target;
    output_struct.t_m = signal.t_m;

    %% ISAR Imaging
    
    % formulate the raw isar image
    [rx_signal, output_struct] = ...
        form_raw_image( signal, target, output_struct);
    p = plotting();
    p.plot_raw_isar(rx_signal, signal.range_array, output_struct.ranges)

    fprintf('Performing post-processing.\n')
    
    % perform range compression via a matched filter by
    % evaluating a convolution of the received signal with a
    % time reverse conjugate of the transmitted signal
    fprintf('   Range compression.\n')
    [rx_signal_range_compressed,output_struct] = ...
        range_compression(...
        signal, ...
        rx_signal, ...
        output_struct);
    
    % perform backprojection to form ISAR image in the target
    % reference frame
    fprintf('   Back-projection.\n')
    [rx_signal_bp, output_struct] = backprojection(...
        signal, ...
        rx_signal_range_compressed, ...
        output_struct.ranges, ...
        output_struct.range_array, ...
        target, ... % assumes height is constant
        yaw, ...
        output_struct.target_positions, ...
        output_struct);
    
    %% Plot results
    fprintf('Simulation complete: %4.2f seconds\n', toc(simStart))
    fprintf('Plotting results.\n')
    
end

