function output_struct = isar_imager(signal,target, sim_params)
    
    %% Preprocssing

    const = Constants();

    simStart = tic;
    
    % save data
    output_struct.signal = signal;
    output_struct.target = target;
    output_struct.t_m = signal.t_m;

    % check input simulation parameters
    if ~sim_params.range_doppler ...
            && ~sim_params.backprojection
        error('Either Range-Doppler processing or Backprojection must be enabled.')
    end

    if sim_params.backprojection
        % calculate the angular extent over the course of the target
        % flyout
        yaw = target.yaw .* signal.t_m ...   
            + (1/2) * target.yawing_rate .* signal.t_m .^ 2;
        angular_extent = yaw(end);
        if angular_extent > 2*pi; angular_extent = 2*pi; end
        signal.cross_range_resolution = ...
            const.c / (2 * angular_extent * signal.fc);
    
        % signal.cross_range_resolution = 1;
        % yaw = zeros(signal.num_pulses,1);
    end

    %% Range Compression
    
    % formulate the raw isar image
    [rx_signal, output_struct] = ...
        form_raw_image( signal, target, output_struct );
    
    fprintf('Performing post-processing.\n')
    
    if sim_params.range_compression

        % perform range compression via a matched filter by
        % evaluating a convolution of the received signal with a
        % time reverse conjugate of the transmitted signal
        fprintf('   Range compression.\n')
        [rx_signal,output_struct] = ...
            range_compression(...
            signal, ...
            rx_signal, ...
            output_struct);
    end
    
    %% Range Alignment

    if sim_params.range_alignment

        % align the range profiles across all pulses relative to a
        % reference range profile
        fprintf('   Range alignment.\n')
        [rx_signal,output_struct] = range_tracking( ...
            signal, ...
            rx_signal, ...
            output_struct);
    end

    %% Phase Adjustment

    if sim_params.phase_adjustment

        % perform coarse phase correction based on dominant
        % scatterers
        fprintf('   Phase adjustment.\n')
        [rx_signal,output_struct] = phase_adjustment(...
            rx_signal, ...
            3, ... % number of bright scatterers
            output_struct);
    end

    %% Image formation

    if sim_params.range_doppler

        % perform Range-Doppler processing across pulses to
        % formulate the final ISAR image
        fprintf('   Range-Doppler.\n')
        [rx_signal,output_struct] = rd_processing( ...
            rx_signal,...
            output_struct);
    
    else
        
        % perform backprojection to form ISAR image in the target
        % reference frame
        fprintf('   Back-projection.\n')
        [rx_signal, output_struct] = backprojection(...
            signal, ...
            rx_signal, ...
            output_struct.ranges, ...
            output_struct.range_array, ...
            target, ... % assumes height is constant
            yaw, ...
            output_struct.target_positions, ...
            output_struct);

    end

    %% Plot results
    fprintf('Simulation complete: %4.2f seconds\n', toc(simStart))
    fprintf('Plotting results.\n')
    
end

