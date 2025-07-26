function [rx_signal_aligned,output_struct] = range_tracking( ...
    signal, ...
    rx_signal_range_compressed, ...
    output_struct)
    
    % Range tracking or range alignment functions to compensate
    % for the translational motion, and aligning the range
    % profiles. A reference range profile is selected along a
    % pulse, and the cross correlation is taken for each pulse.
    % The resulting delay between the reference and each pulse
    % is the amount the range profile needs to be shifted in
    % order to be aligned
    
    % select a reference pulse, the center cross-range bin
    ref_pulse = round(signal.num_pulses / 2);
    
    % extract the corresponding range profile
    ref_profile = rx_signal_range_compressed(ref_pulse, :);
    
    % iterate over each pulse and align the range profile with
    % the reference range profile
    rx_signal_aligned = zeros(size(rx_signal_range_compressed));
    for ipulse = 1:signal.num_pulses
    
        % extract the range profile
        profile = rx_signal_range_compressed(ipulse, :);
    
        % evaluate the cross correlation
        [corr,lag] = xcorr(profile, ref_profile);
    
        % find the lag that maximizes the absolute value of the
        % correlation
        [~, max_idx] = max(abs(corr));
    
        % shift the range profile to maximize correlation with
        % the reference profile
        rx_signal_aligned(ipulse, :) = ...
            circshift(profile, -lag(max_idx));
    
    end
    
    % save data
    output_struct.rx_signal_aligned = rx_signal_aligned;
end

