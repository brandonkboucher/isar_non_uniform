classdef LFM_Signal
    %LFM_SIGNAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fc
        B
        prf
        fs
        Tp

        dt_slow
        dt_fast_time

        t_hat
        t_chirp
        range_array

        range_resolution
        num_pulses
        num_range_bins
        num_chirp_bins


        tx_signal
    end
    
    methods
        function obj = LFM_Signal(fc,B,prf,fs,Tp,T)
            
            % instantiate constants
            const = Constants;
            
            % set radar parameters
            obj.fc = fc; % [Hz]
            obj.B = B; % [Hz]
            obj.prf = prf; % [Hz] this will become non-linear
            obj.fs = fs; % [Hz]
            obj.Tp = Tp; % [s]

            if obj.B * 2 > obj.fs
                warning('Sampling rate does not meet Nyquist criterion.')
            end
            
            % define the chirping rate
            mu = obj.B / obj.Tp;
            
            % slow time, pulse response interval (PRI)
            obj.dt_slow = 1/obj.prf; % [s]
            
            % create a time array corresponding to the fast time
            % sampling
            max_expected_range = 520; % [m]
            tau_max = 2 * max_expected_range / const.c;
            padding = 1 * const.us2s; % [s] additional padding
            
            % define the fast time array
            obj.dt_fast_time = 1/obj.fs;
            obj.t_hat = (0:obj.dt_fast_time:(obj.Tp + tau_max + padding)-(1/obj.fs))';
            
            % define the chirp time array
            obj.t_chirp = (0:obj.dt_fast_time:obj.Tp - (1/obj.fs))';
            obj.range_array = obj.t_hat .* const.c / 2;
            
            % calculate the range resolution
            obj.range_resolution = const.c / (2 * obj.B);
            
            % calculate the number of pulses
            obj.num_pulses = round(T / obj.dt_slow);
            obj.num_range_bins = size(obj.t_hat,1);
            obj.num_chirp_bins = size(obj.t_chirp, 1);

            rect_tx = abs(obj.t_chirp/obj.Tp) <= 1/2;
            obj.tx_signal = ...
                rect_tx .* exp(1j * pi * (mu * obj.t_chirp.^2)); % baseband
            % tx_signal = exp( 2 * pi * -1j * fc * t_hat); % continuous wave

        end
        
    end
end

