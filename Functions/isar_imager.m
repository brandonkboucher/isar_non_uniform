function output_struct = isar_imager(signal,target)
    
    % instantiate constants
    const = Constants;

    simStart = tic;
    
    % save data
    output_struct.signal = signal;
    output_struct.target = target;
    output_struct.t_m = signal.t_m;


    %% ISAR Imaging
    
    % formulate the raw isar image
    [rx_signal, output_struct] = ...
        form_raw_image( signal, target, output_struct );
    
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
    
    % align the range profiles across all pulses relative to a
    % reference range profile
    fprintf('   Range alignment.\n')
    [rx_signal_aligned,output_struct] = range_tracking( ...
        signal, ...
        rx_signal_range_compressed, ...
        output_struct);

    % perform coarse phase correction based on dominant
    % scatterers
    % fprintf('   Phase adjustment.\n')
    % [rx_autofocused,output_struct] = phase_adjustment(...
    %     rx_signal_aligned, ...
    %     3, ... % number of bright scatterers
    %     output_struct);

    % perform Range-Doppler processing across pulses to
    % formulate the final ISAR image
    fprintf('   Range-Doppler.\n')
    [rx_signal_rd,output_struct] = rd_processing( ...
        rx_signal_aligned,...
        output_struct);
    
    %% Plot results
    fprintf('Simulation complete: %4.2f seconds\n', toc(simStart))
    fprintf('Plotting results.\n')
    
end

