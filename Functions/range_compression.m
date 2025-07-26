function [rx_signal_range_compressed,output_struct] = ...
    range_compression(...
    signal, ...
    rx_signal, ...
    output_struct)
    
    const = Constants();

    % range compression via match filtering using the
    % transmitted signal pulse
    h = conj(flipud(signal.tx_signal));
    rx_signal_range_compressed = zeros(size(rx_signal));
    for ipulse =1:signal.num_pulses    
        rx_signal_range_compressed(ipulse,:) = ...
            conv(rx_signal(ipulse,:), h, "same");
    end
    
    % adjust the range axis by the match filter's group delay
    range_array = signal.range_array ...
        - (signal.Tp * const.c / 4);
    output_struct.range_array = range_array;
    
    % save data 
    output_struct.rx_signal_range_compressed = ...
        rx_signal_range_compressed;
end

