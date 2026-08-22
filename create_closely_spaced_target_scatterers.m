function [target_locations, amb_of_k] = create_closely_spaced_target_scatterers( ...
    Ks, ...                 % number of target scatterers
    amb_index, ...          % index of ambiguities (e.g. if 3 -> [-1, 0, 1])
    Wx, ...                 % [m] width of ambiguity 
    Wy, ...                 % [m] array of range values
    separation ...          % [m] separation between target scatterers in crossrange
    )

    % the goal of this script is to create an asymmetric
    % target with close target scatterers. i think there
    % will be some implicit limitation that we should
    % probably abide by like the raylieh resolution, but i
    % am not entirely sure how to calculate that yet. for
    % the time being i think it's fine to just throw
    % something at the wall and see if it sticks.

    num_of_amb = numel(amb_index);

    % ensure target is within bounds
    Wx_eff = Wx - ceil(separation * (Ks+1) / 2);
    Wy_eff = Wy - ceil(separation * (Ks+1) / 2);

    % randomly select a crossrange center of the target for 
    % each ambiguity
    cross_range_center  = amb_index * Wx ...
        + (rand(1,num_of_amb) - 0.5) * Wx_eff;
    
    % randomly pick range center within range
    range_center = (rand(1,num_of_amb) - 0.5) * Wy_eff;
    
    % randomly define an orientation of the target
    theta = 2*pi*rand(num_of_amb,1);
    
    target_locations = zeros(Ks,2);
    
    % round-robin so every requested ambiguity holds at least one scatterer
    amb_of_k = amb_index(mod(0:Ks-1, num_of_amb) + 1).';
    amb_of_k = sort(amb_of_k);

    k = 1;
    for iamb = 1:num_of_amb

        num_scatterers_in_amb = ...
            numel(find(amb_of_k == amb_index(iamb)));

        % target_locations = [cross_range_center + ]
        if mod(num_scatterers_in_amb,2) == 0
            target_idx = -floor(num_scatterers_in_amb/2)...
                :floor(num_scatterers_in_amb/2)-1;
        else
            target_idx = -floor(num_scatterers_in_amb/2)...
                :floor(num_scatterers_in_amb/2);
        end
    
        % position target scatterers along orientation angle,
        % centered at (range_center, cross_range_center)
        cross_range_position = cross_range_center(iamb) ...
            + target_idx * separation * sin(theta(iamb));
        range_position = range_center(iamb) ...
            + target_idx * separation * cos(theta(iamb));
        target_locations(k:num_scatterers_in_amb+k-1, :) = ...
            [cross_range_position.', range_position.'];
        k = k + num_scatterers_in_amb;
    end
end















