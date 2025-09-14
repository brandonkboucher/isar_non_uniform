

classdef Target_L < handle
    %TARGET Summary of this class goes here
    %   Detailed explanation goes here
    
    properties

        % define the state of the center of the target
        position = zeros(1,3) % [m] [x, y, z]
        velocity = zeros(1,3) % [m/s] [x, y, z]
        acceleration = zeros(1,3) % [m/s/s] [x, y, z]
        velocity_magnitude = 0; % [m/s]
    
        % if the target is traveling along a circular path
        radius = 0; % [m]
        yawing_rate = 0; % [rad/s]
        period = 0; % [s]
        dyaw = 0; % [rad]
        circle_center = zeros(1,2);
        scatterer_angle_to_heading = [];

        % define orientation
        yaw = 0; % [rad]
        yaws = []; % [rad] array of yaws

        % define the number of scatters and the position of
        % each scatter
        num_scatterers = 3;
        scatter_positions
        scatter_velocities
        scatter_configuration
        scatter_relative_positions

        dt = 0.01; % [s] time iteration
    end
    
    methods

        function obj = Target_L( ...
                dt, ...
                initial_center_position, ...
                yaw)
            
            % the previous basic target model allowed for
            % the propagation of a target with three
            % scatters in a straight line. We now want the
            % target to rotate along the y axis with angular
            % acceleration. many of the terms between the
            % two models would overlap, so I'm going to
            % separate them entirely to minimize confusion

            obj.dt = dt;
            obj.position = initial_center_position;

            % currently the target configuration has three
            % scatters that are initialized as:

            % y-axis
            % ^
            % |               x (point 1)
            % |               ^
            % | target length |
            % |               |
            % |               |
            % |               |
            % |               |     r = radius to scatter
            % |     x         |----------->x
            % | (point 2)             (point 3)
            % |---------------------------------> x-axis

            % where the target is translationally moving
            % along the y-axis

            % when rolling in the counterclock direction,
            % point 2 is moving positively such that the
            % z-component is increasing

            % define the shape of the simulated aircraft
            % as a triangle relative to the center
            obj.radius = [3, 7];
            
            % set the initial yaw and define a rotation
            % matrix
            obj.yaw = yaw;
            obj.yaws = [obj.yaws; obj.yaw];
            R_yaw = ...
            [cos(obj.yaw), -sin(obj.yaw), 0; ...
             sin(obj.yaw), cos(obj.yaw), 0; ...
             0, 0, 1];

            obj.scatter_configuration = ...
                [0, 0, 0; ... % scatter 1
                obj.radius(2), 0, 0; ... % scatter 2
                0, obj.radius(1), 0];     % scatter 3

            % the scatterer's positions relative to the
            % initial scatterer configuration
            obj.scatter_relative_positions = ...
                obj.scatter_configuration;

            % transform to the radar's reference frame
            obj.scatter_positions = ...
                (R_yaw * obj.scatter_configuration')' ...
                + obj.position;

        end
        
        function obj = propagate(obj)

            % update the roll using the body-axis roll rate
            % and body-axis acceleration, please note that
            % this is a body axis angle
            obj.yaw = obj.yaw + obj.yawing_rate * obj.dt;
            obj.yaws = [obj.yaws; obj.yaw];

            % define a rotation matrix about the z-axis
            rot_z = ...
            [cos(obj.yaw), -sin(obj.yaw),0; ... 
            sin(obj.yaw), cos(obj.yaw), 0; ... 
            0, 0, 1];

            % rotate the scatterers relative to the center
            % of the target
            obj.scatter_relative_positions = ...
                (rot_z * (obj.scatter_configuration'))';

            obj.scatter_positions = ...
                obj.scatter_relative_positions ...
                + obj.position;

            obj.set_scattering_velocities();
        end

        function set.yawing_rate(obj,yawing_rate)
            obj.yawing_rate = yawing_rate;
            obj.set_scattering_velocities();
           
        end

        function set_scattering_velocities(obj)

            for iscatterer = 1:obj.num_scatterers
                obj.scatter_velocities(iscatterer, :) = cross( ...
                    obj.scatter_relative_positions(iscatterer, :), ...
                    [0,0,obj.yawing_rate]);
            end
        end

    end
end

