

classdef Target_Simple < handle
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

        euler_angles = zeros(1,3); % [rad] roll, pitch, yaw
        body_axis_rotation_vel = zeros(1,3); % [rad/s] p, q, r
        body_axis_rotation_acc = zeros(1,3); % [rad/s/s] p dot, q dot, r dot
      
        % define the number of scatters and the position of
        % each scatter
        num_scatterers = 3;
        scatter_positions
        scatter_velocities
        scatter_configuration
        scatter_relative_positions

        dt; % [s] time iteration
    end
    
    methods

        function obj = Target_Simple(dt, initial_position, num_scatterers)
            
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

            % if obj.num_scatterers > 1
            %     warning('The Target_Simple model can only handle one target scattering point')
            %     obj.num_scatterers = 1;
            % end

            obj.radius_to_scatterer = 0;
            obj.scatter_configuration = [0, obj.radius_to_scatterer, 0]; 
            
            % for now we will assume that for the simple
            % target that the target is placed at x = z = 0,
            if num_scatterers > 1
                if mod(num_scatterers, 2) == 0
                    line_separation = 10;
    
                    yline_array =  ...
                    - line_separation * (num_scatterers - 1): ...
                    line_separation:line_separation * (num_scatterers - 1);
    
                else
                    line_separation = 10;
    
                    yline_array =  ...
                    - line_separation * (num_scatterers - 2): ...
                    line_separation:line_separation * (num_scatterers - 2);
    
                end
    
                obj.scatter_configuration = zeros(num_scatterers, 3);
                obj.scatter_configuration(:,2) = yline_array';
            else
                obj.scatter_configuration = zeros(1,3);
            end
            % the scatterer's positions relative to the
            % initial scatterer configuration
            obj.scatter_relative_positions = ...
                obj.scatter_configuration;

            % transform to the radar's reference frame
            obj.scatter_positions = ...
                obj.scatter_configuration ...
                + repmat(obj.position, num_scatterers, 1);

            % there is no translation or tangential velocity
            % as the radius of the scatter is zero (point 
            % target)
            obj.scatter_velocities = zeros(num_scatterers, 3);

        end
        
        function obj = propagate(obj)

            obj.update_rotational_dynamics()
            obj.update_translational_dynamics()

        end

        function update_rotational_dynamics(obj)

            % update the body-axis rotational rate
            obj.body_axis_rotation_vel = ...
                obj.body_axis_rotation_vel ...
                + obj.body_axis_rotation_acc * obj.dt;

            % update the euler angles
            obj.euler_angles = obj.euler_angles ...
                + obj.body_axis_rotation_vel * obj.dt ...
                + (1/2) * obj.body_axis_rotation_acc * obj.dt^2;

            % define the rotation matrix about the z-axis
            %rot_z = [cos()]

            % the roll of the target modifies the position
            % of the back two scatterers in addition to the
            % translational motion


        end

        function update_translational_dynamics(obj)

            % propagate target velocity
            obj.velocity = obj.velocity ...
                + obj.acceleration * obj.dt; 
            obj.scatter_velocities = ...
                repmat(obj.velocity, obj.num_scatterers, 1);

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

