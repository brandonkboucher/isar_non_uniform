classdef scenario_basic < handle
    

    
    properties
        signal
        target
    end
    
    methods
        function obj = scenario_basic(T)
            
            % instantiate constants
            const = Constants;
            c = const.c;

            % simulation duration
            T = 1; % [s]
                        
            % initialize radar parameters and LFM signal model
            fc  = 10 * const.GHz2Hz; % [Hz] center frequency - X-band
            B   = 149.9 * const.MHz2Hz; % [Hz] bandwidth
            prf = 1 * const.kHz2Hz; % [Hz] pulse repetition frequency
            fs  = 600 * const.MHz2Hz; % [Hz] sampling frequency
            Tp  = 10 * const.us2s; % [s] pulse width

            % define the LFM signal
            obj.signal = Signal_Basic(fc, B, prf, fs, Tp, T);

            % define rotation rate
            yawing_rate = pi/4;

            % intialize target
            obj.target = Target_Basic(...
                obj.signal.t_slow, ...
                yawing_rate);
        end
       
    end
end

