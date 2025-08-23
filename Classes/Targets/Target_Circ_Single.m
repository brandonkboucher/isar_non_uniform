

classdef Target_Circ_Single < handle
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

        % define the number of scatters and the position of
        % each scatter
        num_scatterers = 1;
        scatter_positions
        scatter_velocities
        scatter_configuration
        scatter_relative_positions

        dt = 0.01; % [s] time iteration
    end
    
    methods

        function obj = Target_Circ_Single(dt, initial_position, num_scatterers)
            
            % the previous basic target model allowed for
            % the propagation of a target with three
            % scatters in a straight line. We now want the
            % target to rotate along the y axis with angular
            % acceleration. many of the terms between the
            % two models would overlap, so I'm going to
            % separate them entirely to minimize confusion

            obj.dt = dt;
            obj.position = initial_position;
            obj.num_scatterers = num_scatterers;

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
            obj.radius = 20;
            
            obj.scatter_configuration = ...
                [obj.radius,0,0];                 % scatter 3

            % the scatterer's positions relative to the
            % initial scatterer configuration
            obj.scatter_relative_positions = ...
                obj.scatter_configuration;

            % transform to the radar's reference frame
            obj.scatter_positions = ...
                obj.scatter_configuration + obj.position;

        end
        
        function obj = propagate(obj)

            % update the roll using the body-axis roll rate
            % and body-axis acceleration, please note that
            % this is a body axis angle
            obj.yaw = obj.yaw + obj.yawing_rate * obj.dt;

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

        % function set.velocity_magnitude(obj, velocity_magnitude)
        %     obj.velocity_magnitude = velocity_magnitude;
        % 
        % end
        % 
        % function obj = calculate_tangential_velocity(obj)
        % 
        %     % using the position vector and the rotation
        %     % rate, calculate the tangential velocity
        %     r_xy = obj.position(1:2) - obj.circle_center;
        %     obj.velocity = cross(...
        %         [0,0, obj.yawing_rate], [r_xy, 0]);
        %     obj.velocity_magnitude = norm(obj.velocity);
        % end
        % 
        % function set.radius(obj, radius)
        % 
        %     obj.radius = radius;
        % 
        %     if obj.yawing_rate == 0
        %         % if the radius is set, the trajectory is
        %         % assumed circular, set the related properties
        %         obj.initialize_circular_trajectory_from_velociy();
        %     else
        %         obj.initialize_circular_trajectory_from_omega();
        %     end
        % 
        % end
        % 
        % function obj = initialize_circular_trajectory_from_velociy(obj)
        % 
        %     if norm(obj.velocity_magnitude) <= 0
        % 
        %         warning('Please set the velocity prior to defining the rotation parameters.')
        %     else
        %         % calculate the rotation rate
        %         obj.yawing_rate = obj.velocity_magnitude / obj.radius; % [rad/s]
        %         obj.period = 2 * pi / obj.yawing_rate; % [s]
        % 
        %         % for each time interval, calculate the 2D angle
        %         % change along the circular path
        %         obj.dyaw = obj.yawing_rate * obj.dt;
        % 
        %         % set the circle's center points
        %         obj.circle_center = ...
        %             [obj.position(1) - obj.radius * cos(obj.yaw), ...
        %             obj.position(2) - obj.radius * sin(obj.yaw)];
        % 
        %         % calculate tangential velocity
        %         obj.calculate_tangential_velocity();
        % 
        %     end
        % 
        % end

        % function obj = initialize_circular_trajectory_from_omega(obj)
        % 
        %     % calculate the rotation rate
        %     % obj.omega = obj.velocity_magnitude / obj.radius; % [rad/s]
        %     obj.period = 2 * pi / obj.yawing_rate; % [s]
        % 
        %     % for each time interval, calculate the 2D angle
        %     % change along the circular path
        %     obj.dyaw = obj.yawing_rate * obj.dt;
        % 
        %     % set the circle's center points
        %     obj.circle_center = ...
        %         [obj.position(1) - obj.radius * cos(obj.yaw), ...
        %         obj.position(2) - obj.radius * sin(obj.yaw)];
        % 
        %     % calculate tangential velocity
        %     obj.calculate_tangential_velocity();
        % 
        % 
        % end

    end
end

