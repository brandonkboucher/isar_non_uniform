classdef scenario_complex_traj_jerk < handle
    
    properties
        signal
        target

        angular_extent = nan
        T = nan
    end
    
    methods
        function obj = scenario_complex_traj_jerk(T, angular_extent)
            
            obj.T = T;
            obj.angular_extent = angular_extent;

            if ~isnan(T) && ~isnan(angular_extent)
                error('The simulation duration and angular extent of the target cannot be simultaneously set.')
            end

            const = Constants();

            % initialize radar parameters and LFM signal model
            fc  = 1 * const.GHz2Hz; % [Hz] center frequency - X-band
            B   = 149.9 * const.MHz2Hz; % [Hz] bandwidth
            prf = 10 * const.kHz2Hz; % [Hz] pulse repetition frequency
            fs  = 600 * const.MHz2Hz; % [Hz] sampling frequency
            Tp  = 10 * const.us2s; % [s] pulse width

            % define the LFM signal
            obj.signal = Signal_Basic(fc, B, prf, fs, Tp, obj.T);

            % intialize target
            obj.target = Target_complex_traj_jerk(obj.signal.t_slow);

            fprintf('Simulation duration:       %2.2f seconds \n', obj.T)
            fprintf('Angular extent:            %2.2f degress\n',...
                rad2deg(max(obj.target.yaws)))
            % fprintf('Angular velocity:          %2.2f rad/s\n', yawing_rate)

        end

        function check_parameters(obj)
        
            % check Nyquist criterion
            if obj.signal.B * 2 > obj.signal.fs
                warning('Sampling rate does not meet Nyquist criterion.')
            end

            % instantiate constants
            const = Constants;
            c = const.c;

            % check the crossrange is adequately sampled
            lambda = c / obj.signal.fc;
            angular_extent = obj.yaws(end);
            
            % calculate the length of the aperture
            L = obj.target.target_center_position ...
                * angular_extent;

            % calculate the angular sampling
            delta_theta = lambda / (2 * L);

            % calculate the number of samples per aperture
            %N_samples = angular_extent / delta_theta;

            % calculate the unambiguous range
            Ramb = c / (2 * obj.signal.prf);

            if obj.target.target_center_position > Ramb
                warning('The target range is ambiguous.')
            end

            if obj.signal.Tp * obj.signal.prf > 0.05
                warning('Duty cycle is too high')
            end

        end

        function obj = reset_duration(obj)
            
            % obj.T = calculate_simulation_duration(...
            %         obj.angular_extent, ...
            %         obj.target.yawing_rate, ...
            %         obj.target.yawing_acceleration);

            obj.T = obj.angular_extent / (pi/4);

            yawing_rate = calculate_variable_omega(...
                pi/4, ...
                obj.target.yawing_acceleration, ...
                obj.T);

            % define the LFM signal
            obj.signal = Signal_Basic( ...
                obj.signal.fc, ...
                obj.signal.B, ...
                obj.signal.prf, ...
                obj.signal.fs, ...
                obj.signal.Tp, ...
                obj.T);

            % intialize target
            obj.target = Target_Basic(...
                obj.signal.t_slow, ...
                yawing_rate, ...
                obj.target.yawing_acceleration);

            fprintf('Simulation duration:       %2.2f seconds \n', obj.T)
            fprintf('Angular extent:            %2.2f degress\n',...
                rad2deg(obj.target.yaws(end)))
            fprintf('Angular velocity:          %2.2f rad/s\n', yawing_rate)
        end
    end
end

