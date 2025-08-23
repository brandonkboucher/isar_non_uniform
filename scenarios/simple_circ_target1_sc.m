classdef simple_circ_target1_sc < handle
    
    % the goal of this scenario is a uniformly,
    % non-accelerating target without translational motion.
    
    
    properties
        signal
        target
    end
    
    methods
        function obj = simple_circ_target1_sc(T)
            
            % instantiate constants
            const = Constants;

            % define the number of scattering points off of the target
            num_scatterers = 1;
            
            % define the target's initial position
            target_position = [0, 100, 0];
            
            % initialize radar parameters and LFM signal model
            fc  = 10 * const.GHz2Hz; % [Hz] center frequency - X-band
            B   = 149.9 * const.MHz2Hz; % [Hz] bandwidth
            prf = 1 * const.kHz2Hz; % [Hz] pulse repetition frequency
            fs  = 300 * const.MHz2Hz; % [Hz] sampling frequency
            Tp  = 5 * const.us2s; % [s] pulse width
            
            % in order to pad the fast time array, a rough guess at the
            % maximum range must be set
            max_range = norm(target_position);
            
            % define the LFM signal
            obj.signal = LFM_Signal(fc, B, prf, fs, Tp, T, max_range);
                
            % intialize target
            obj.target = Target_Circ_Single(...
                obj.signal.dt_slow, target_position, num_scatterers);
            
            % set straight line velocity towards the radar
            obj.target.yawing_rate = pi; 

        end
        
        
    end
end

