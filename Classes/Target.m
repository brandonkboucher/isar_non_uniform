

classdef Target < handle
    %TARGET Summary of this class goes here
    %   Detailed explanation goes here
    
    properties

        % translational dynamics
        position = zeros(1,3) % [m] [x, y, z]
        velocity = zeros(1,3) % [m/s] [x, y, z]
        acceleration = zeros(1,3) % [m/s/s] [x, y, z]
        velocity_magnitude = 0; % [m/s]
    
        % rotational dynamics
        radius_to_scatterer = 0; % [m]
        target_length = 0;
        p = 0; % [rad/s] roll rate
        p_dot = 0; % [rad/s/s] derivative of roll rate
        roll = 0;

        % define the number of scatters and the position of
        % each scatter
        num_scatters = 3;
        scatter_positions
        scatter_configuration
        scatter_relative_positions

        dt; % [s] time iteration
    end
    
    methods

        function obj = Target(dt, initial_position, num_scatters)
            
            % the previous basic target model allowed for
            % the propagation of a target with three
            % scatters in a straight line. We now want the
            % target to rotate along the y axis with angular
            % acceleration. many of the terms between the
            % two models would overlap, so I'm going to
            % separate them entirely to minimize confusion

            obj.dt = dt;
            obj.position = initial_position;
            obj.num_scatters = num_scatters;

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
            obj.radius_to_scatterer = 5;
            obj.target_length = 10;

            obj.scatter_configuration = ...
                [0, -obj.target_length, 0; ... % scatter 1
                -obj.radius_to_scatterer, 0, 0; ... % scatter 2
                obj.radius_to_scatterer, 0, 0]; % scatter 3

            % the scatterer's positions relative to the
            % initial scatterer configuration
            obj.scatter_relative_positions = ...
                obj.scatter_configuration;

            % transform to the radar's reference frame
            obj.scatter_positions = ...
                obj.scatter_configuration + obj.position;

        end
        
        function obj = propagate(obj)

            obj.update_rotational_dynamics()
            obj.update_translational_dynamics()

        end

        function update_rotational_dynamics(obj)

            % update the roll rate
            obj.p = obj.p + obj.p_dot * obj.dt;

            % update the roll using the body-axis roll rate
            % and body-axis acceleration, please note that
            % this is a body axis angle
            obj.roll = obj.roll ...
                + obj.p * obj.dt ...
                + (1/2) * obj.p_dot * (obj.dt)^2;

            % the roll of the target modifies the position
            % of the back two scatterers in addition to the
            % translational motion
            obj.scatter_relative_positions = ...
                [0,-obj.target_length,0; ... % scatterer 1
                -obj.radius_to_scatterer * cos(obj.roll), 0, -obj.radius_to_scatterer * sin(obj.roll); ... % scatterer 2
                obj.radius_to_scatterer * cos(obj.roll), 0, obj.radius_to_scatterer * sin(obj.roll)];

        end

        function update_translational_dynamics(obj)

            % propagate target velocity
            obj.velocity = obj.velocity ...
                + obj.acceleration * obj.dt; 

            % propagate target position
            obj.position = obj.position ...
                + obj.velocity * obj.dt ...
                + (1/2) * obj.acceleration * (obj.dt)^2;

            % transform to the radar's reference frame
            obj.scatter_positions = ...
                obj.scatter_relative_positions ...
                + obj.position;

        end

        function set.velocity_magnitude(obj, velocity_magnitude)
            obj.velocity_magnitude = velocity_magnitude;

        end
   
    end
end

