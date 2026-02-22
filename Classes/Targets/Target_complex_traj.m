classdef Target_complex_traj < handle
    %TARGET_BASIC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        num_scatterers = 3;
        target_center_position = [0, 1000, 0];
        scatter_relative_positions = [7, 0, 0; ...
                                      0, 3, 0; ...
                                      0, 0, 0];
        yawing_rate = pi/4;
        yawing_acceleration = 0;
        yaws

        t_slow

    end
    
    methods
        function obj = Target_complex_traj(t_slow)

            % varargin must be of the form:
            % 1. num_scatterers
            % 2. target_center_position
            % 3. scatter_relative_positions

            obj.t_slow = t_slow;

            % initial motion parameters
            initial_yaw = 0;
            initial_yawing_rate = pi/8;

            yawing_acceleration = 2*pi*[1, -1, -1, 1];
            num_partitions = size(yawing_acceleration,2);

            % partition into four segments
            idx_segment_incr = ...
                round(size(t_slow,1)/num_partitions);
            
            yaws = zeros(size(t_slow,1),1);
            yawing_rate = zeros(size(t_slow,1),1);
            yawing_acc = zeros(size(t_slow,1),1);

            for iseg = 1:num_partitions

                % define the segment indices
                idx = 1 + (iseg-1)*idx_segment_incr...
                    :iseg*idx_segment_incr;

                if iseg == num_partitions
                    if idx(end) > size(t_slow,1)
                        idx = idx(idx <= size(t_slow,1));
                    
                    elseif idx(end) < size(t_slow,1)
                        missing_idx = 1:(size(t_slow,1) - idx(end));
                        idx = [idx, idx(end) + missing_idx];
                    end
                end

                t_seg = t_slow(idx) - t_slow(idx(1));

                yaws(idx) = ...
                    initial_yaw ...
                    + initial_yawing_rate * t_seg ...
                    + (1/2) * yawing_acceleration(iseg) * t_seg .* t_seg;
                
                yawing_rate(idx) = ...
                    initial_yawing_rate ...
                    + yawing_acceleration(iseg) .* t_seg;
                yawing_acc(idx) = repmat(yawing_acceleration(iseg), size(idx,1),1);
                
                initial_yawing_rate = yawing_rate(idx(end));
                initial_yaw = yaws(idx(end));
                
            end
            obj.yaws = yaws;
            obj.yawing_rate = yawing_rate;
            obj.yawing_acceleration = yawing_acc;
        end
    end
end

