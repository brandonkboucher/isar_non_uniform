
function [rx_autofocused,output_struct] = phase_adjustment(...
    rx_signal_aligned, ...
    K, ... % number of bright scatterers
    output_struct)
    
    %% Phase Adjustment
    
    % Phase adjustment is a way to sharpen the image along the
    % cross-range direction prior to Range-Doppler processing,
    % and works by assuming the phase of the dominant scatters
    % is roughly constant. Any phase change across pulses is
    % assumed error or phase drift. To correct, find the
    % dominant scatters phase angles and average the phases for
    % each pulse.
    
    % I'm having a real tough time finding a coherent, fully 
    % written out version of the MVM so I'm going to rely on the
    % AI to help for now. If I can't get this to work then I'll
    % try implementing the phase gradient algorithm
    
    % find the across all pulses the max absolute of the
    % signal for each range bin
    range_bins_max = max(abs(rx_signal_aligned), [], 1);
    
    % find the K range bins with the max absolute of the
    % received signal. You are finding across all pulses the
    % maximum received signal for K range bins.
    [~, dominant_idx] = maxk(range_bins_max, K);
    dominant_ranges = output_struct.range_array(dominant_idx);

    % extract the phases of the most dominant r
    phases = angle(rx_signal_aligned(:,dominant_idx));
    
    % apply the phase correction over the range bins
    correction = exp(-1j * phases);
    rx_autofocused = rx_signal_aligned .* mean(correction, 2);
    
    % save data
    output_struct.rx_autofocused = rx_autofocused;
end

