classdef Simulation_Parameters < handle
    
    % basic simulation parameter object that holds booleans
    % of what post-processing should be performed
    
    properties
        range_compression = true;
        translational_moco = true;
        range_alignment = true;
        phase_adjustment = true;
        rotational_moco = true;
        backprojection = false;
        range_doppler = false;
    end
    
    methods
        function obj = Simulation_Parameters()
            
        end
        
        function obj = set.translational_moco(obj, translational_moco)

            obj.translational_moco = translational_moco;
            
            if ~translational_moco
                obj.range_alignment = false;
                obj.phase_adjustment = false;
            end

        end

        function obj = set.range_doppler(obj, range_doppler)

            obj.range_doppler   = range_doppler;
        end

        function obj = set.backprojection(obj, backprojection)

            obj.backprojection = backprojection;
            obj.translational_moco = false;
        end
    end
end

