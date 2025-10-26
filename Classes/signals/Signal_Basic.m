classdef Signal_Basic
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
        function obj = Signal_Basic( ...
                fc, ...
                B, ...
                prf, ...
                fs, ...
                Tp, ...
                T)
            
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
            dt_slow = 1/obj.prf; % [s]
            obj.t_slow = (0:dt_slow:T-dt_slow)';
            
            % define the fast time array
            dt_fast_time = 1/obj.fs;
            obj.t_fast = (0:dt_fast_time:(obj.Tp)-(1/obj.fs))';
            
            % define the chirp time array
            obj.range_array = obj.t_fast .* const.c / 2;
            
            % calculate the range resolution
            obj.range_resolution = const.c / (2 * obj.B);
            
            % calculate the number of pulses
            obj.num_pulses = round(T / dt_slow);
            obj.num_range_bins = size(obj.t_fast,1);

        end
        
    end
end

