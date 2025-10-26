classdef Target_Basic
    %TARGET_BASIC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        num_scatterers = 3;
        target_center_position = [0, 1000, 0];
        scatter_relative_positions = [7, 0, 0; ...
                                      0, 3, 0; ...
                                      0, 0, 0];
        yawing_rate = pi/4;
        yaws

    end
    
    methods
        function obj = Target_Basic(...
                t_slow,...
                yawing_rate,...
                varargin)

            % varargin must be of the form:
            % 1. num_scatterers
            % 2. target_center_position
            % 3. scatter_relative_positions

            if ~isempty(varargin)

                % number of scatterers
                obj.num_scatterers = varargin{1};
                
                % define the target's initial position
                obj.target_center_position = ...
                    varargin{2};
                
                obj.scatter_relative_positions = ...
                    varargin{3};
            end
            
            obj.yawing_rate = yawing_rate;

            % calculate the yaws
            obj.yaws = obj.yawing_rate * t_slow;
        end
        
    end
end

