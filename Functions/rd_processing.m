function [rx_signal_rd,output_struct] = rd_processing( ...
    rx_autofocused,...
    output_struct)
    
    
    % perform Doppler processing (or az FFT) on the range
    % compressed data   
    rx_signal_rd = ...
        fftshift(fft(rx_autofocused, [], 1), 1);
    
    % save data
    output_struct.rx_signal_rd = rx_signal_rd;
end

