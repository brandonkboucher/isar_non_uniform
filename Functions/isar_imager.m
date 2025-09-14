function output_struct = isar_imager(signal,target, sim_params)
    
    %% Preprocssing

    const = Constants();

    simStart = tic;
    
    % save data
    output_struct.signal = signal;
    output_struct.target = target;
    output_struct.t_slow = signal.t_slow;

    % check input simulation parameters
    if ~sim_params.range_doppler ...
            && ~sim_params.backprojection
        error('Either Range-Doppler processing or Backprojection must be enabled.')
    end

    %% Range Compression
    
    if isa(signal, 'LFM_Signal')

        % formulate the raw isar image using a transmitted
        % LFM signal
        [rx_signal, output_struct] = ...
            form_raw_image_LFM( ...
            signal, ...
            target, ...
            output_struct );
        
    elseif isa(signal, 'Basic_Signal')

        % formulate the raw isar image using a basic signal
        % model approximation
        [rx_signal, output_struct] = ...
            form_raw_image_basic( ...
            signal, ...
            target, ...
            output_struct );
    else
        error('Signal model not set correctly.')
    end

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
        
        % calculate the angular extent of the taget, if
        % greater than 360 degrees, set to 360 degrees
        angular_extent = target.yaws(end);
        if angular_extent > 2*pi; angular_extent = 2*pi; end
        
        % define the cross range resolution for the grid
        signal.cross_range_resolution = ...
            const.c / (2 * angular_extent * signal.fc);
        output_struct.signal.cross_range_resolution =...
            signal.cross_range_resolution;

        % backprojection grid extent
        grid_extent = [20, 20]; % [cross range, range]

        % perform backprojection to form ISAR image in the target
        % reference frame
        fprintf('   Back-projection.\n')
        [rx_signal, output_struct] = backprojection(...
            signal, ...
            rx_signal, ...
            output_struct.ranges, ...
            output_struct.range_array, ...
            target, ... % assumes height is constant
            target.yaws, ...
            output_struct.target_positions, ...
            output_struct, ...
            grid_extent);

        plot_bp_traj(...
            grid_extent, ...
            target.yaws, ...
            output_struct)

    end

    %% Plot results
    fprintf('Simulation complete: %4.2f seconds\n', toc(simStart))
    fprintf('Plotting results.\n')
    
end

