classdef Basic_Signal
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

        t_fast
        t_slow
        t_chirp
        range_array
        doppler_array

        cross_range_resolution
        range_resolution
        num_pulses
        num_range_bins
        num_chirp_bins


        tx_signal
    end
    
    methods
        function obj = Basic_Signal(fc,B,prf,fs,Tp,T, max_expected_range)
            
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
            
            % slow time, pulse response interval (PRI)
            obj.dt_slow = 1/obj.prf; % [s]
            
            % create a time array corresponding to the fast time
            % sampling
            tau_max = 2 * max_expected_range / const.c;
            padding = 1 * const.us2s; % [s] additional padding
            
            % define the fast time array
            obj.dt_fast_time = 1/obj.fs;
            obj.t_fast = (0:obj.dt_fast_time:(obj.Tp)-(1/obj.fs))';

            % slow time array
            obj.t_slow = (0:obj.dt_slow:T-obj.dt_slow)';

            % define the chirp time array
            obj.t_chirp = (0:obj.dt_fast_time:obj.Tp - (1/obj.fs))';
            obj.range_array = obj.t_fast .* const.c / 2;
            obj.doppler_array = linspace(-obj.prf/2, obj.prf/2, size(obj.t_slow,1))';

            % calculate the range resolution
            obj.range_resolution = const.c / (2 * obj.B);
            
            % calculate the number of pulses
            obj.num_pulses = round(T / obj.dt_slow);
            obj.num_range_bins = size(obj.t_fast,1);
            obj.num_chirp_bins = size(obj.t_chirp, 1);

            t_tx = (-obj.Tp/2 : obj.dt_fast_time : obj.Tp/2 - obj.dt_fast_time)';
            obj.tx_signal = sinc(obj.B * t_tx); % baseband

        end
        
    end
end

