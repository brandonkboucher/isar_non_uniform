

classdef Target < handle
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
        omega = 0; % [rad/s]
        period = 0; % [s]
        dtheta = 0; % [rad]
        circle_center = zeros(1,2);

        % define orientation
        heading = 0; % [rad]

        % define the number of scatters and the position of
        % each scatter
        num_scatters = 1;
        scatter_positions
        scatter_configuration

        dt = 0.01; % [s] time iteration
    end
    
    methods

        function obj = Target(dt, initial_position, num_scatters)
            
            obj.dt = dt;
            obj.position = initial_position;
            obj.num_scatters = num_scatters;

            % using the number of scatter points define the
            % location of each point relative to the center
            % of the target
            
            if obj.num_scatters == 1

                % if there is only one scatter point then
                % place it at the center of the target
                obj.scatter_configuration = [0,0,0];
                obj.scatter_positions = obj.position;
            
            elseif obj.num_scatters == 3
                
                % define the shape of the simulated aircraft
                % as a triangle relative to the center
                obj.scatter_configuration = ...
                    [0, 20, 0; ... % scatter 1
                    10, -20, 0; ... % scatter 2
                    -10, -20, 0]; % scatter 3

                % transform to the radar's reference frame
                obj.scatter_positions = ...
                    obj.scatter_configuration + obj.position;

            else
                warning('Currently this simulation is only capable of modeling 1 or 3 scatter points.')
            end

            obj.heading = atan(obj.position(1)/obj.position(2));

        end
        
        function obj = propagate(obj)

            if obj.radius > 0 % circular path

                % change heading using angle change due to
                % yawing
                obj.heading = obj.heading + obj.dtheta; % [rad]
    
                % calculate new target position along circular
                % path
                obj.position(1) = ...
                    obj.circle_center(1) ...
                    + obj.radius * cos(obj.heading);
                obj.position(2) = ...
                    obj.circle_center(2) ...
                    + obj.radius * sin(obj.heading);
                obj.scatter_positions = obj.position;

                % calculate velocity vector
                obj.calculate_tangential_velocity();
            
            else % straight line

                % propagate target position
                obj.position = obj.position ...
                    + obj.velocity * obj.dt ...
                    + (1/2) * obj.acceleration * (obj.dt)^2;
    
                % propagate target velocity
                obj.velocity = obj.velocity ...
                    + obj.acceleration * obj.dt; 

                % transform to the radar's reference frame
                obj.scatter_positions = ...
                    obj.scatter_configuration + obj.position;
            end
        end

        function set.velocity_magnitude(obj, velocity_magnitude)
            obj.velocity_magnitude = velocity_magnitude;

        end

        function obj = calculate_tangential_velocity(obj)

            % using the position vector and the rotation
            % rate, calculate the tangential velocity
            r_xy = obj.position(1:2) - obj.circle_center;
            obj.velocity = cross(...
                [0,0, obj.omega], [r_xy, 0]);
            obj.velocity_magnitude = norm(obj.velocity);
        end

        function set.radius(obj, radius)

            obj.radius = radius;

            % if the radius is set, the trajectory is
            % assumed circular, set the related properties
            obj.initialize_circular_trajectory();

        end

        function obj = initialize_circular_trajectory(obj)

            if norm(obj.velocity_magnitude) <= 0

                warning('Please set the velocity prior to defining the rotation parameters.')
            else
                % calculate the rotation rate
                obj.omega = obj.velocity_magnitude / obj.radius; % [rad/s]
                obj.period = 2 * pi / obj.omega; % [s]

                % for each time interval, calculate the 2D angle
                % change along the circular path
                obj.dtheta = obj.omega * obj.dt;

                % set the circle's center points
                obj.circle_center = ...
                    [obj.position(1) - obj.radius * cos(obj.heading), ...
                    obj.position(2) - obj.radius * sin(obj.heading)];

                % calculate tangential velocity
                obj.calculate_tangential_velocity();

            end

        end

    end
end

