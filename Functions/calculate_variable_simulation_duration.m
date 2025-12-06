function T = calculate_simulation_duration(...
    angular_extent,...% [rad]
    yawing_rate, ... % [rad/s]
    yawing_acceleration ... % [rad/s/s]
    )
    
    % formulate the quadratic equation for the rotational
    % motion of the target
    syms t
    eqn = angular_extent - yawing_rate * t - (1/2)*yawing_acceleration * t^2 == 0;
    
    % solve for the simulation duration required
    S = solve(eqn);
    S = eval(S);

    % choose the positive T
    if size(S,1) > 1 && all(S > 0)
        warning('Both solutions to the quadratic simulation duration equation are positive. Selecting the largest value.')
    end
    T = max(S(S>0));
end

