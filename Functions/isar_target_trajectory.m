function traj = isar_target_trajectory(cfg, opts)
% ISAR_TARGET_TRAJECTORY  Slow-time sampling and target yaw for a scenario.
%
%   traj = ISAR_TARGET_TRAJECTORY(cfg, opts) returns the yaw angle theta(t)
%   the target sweeps out during the collection, together with the slow-time
%   grid it is sampled on and the angular rate that produced it.
%
%   Split out of ISAR_RUN_SCENARIO so the GUI can preview the trajectory
%   without paying for a reconstruction, and so the previewed motion is by
%   construction the motion the simulation uses.
%
%   traj fields:
%     .t_m         [M x 1] slow-time samples [s]
%     .theta       [M x 1] yaw angle [rad]
%     .theta_rate  [M x 1] yaw rate [rad/s], exact for both trajectory models
%     .w0 .w1 .w2  the rate/acceleration/jerk actually in force (w1 and w2
%                  are zeroed unless the scenario is Accelerating)
%     .T           [s] collection duration
%     .maneuvering true when the piecewise maneuvering trajectory was used

    if nargin < 2 || isempty(opts)
        opts = struct();
    end
    opts = isar_default_options(opts);

    T   = (1/opts.prf) * opts.Nd;                        % [s] duration
    t_m = (0:(1/opts.prf):(T - 1/opts.prf)).';           % [s] slow-time

    is_accel = strcmpi(cfg.AngleRate, 'Accelerating');

    w0 = opts.w0;   % [rad/s] yawing rate
    if is_accel
        w1 = opts.w1;   % [rad/s/s] yawing acceleration
        w2 = opts.w2;   % [rad/s/s/s] yawing jerk
    else
        w1 = 0;
        w2 = 0;
    end

    maneuvering = opts.complex_maneuver && is_accel;

    if maneuvering
        [theta, theta_rate] = create_complex_target_trajectory( ...
            w0, opts.w1, opts.w2, opts.jerk_mag, t_m);
    else
        theta = w0 * t_m ...
            + (1/2) * w1 * t_m .* t_m ...
            + (1/3) * w2 * t_m .* t_m .* t_m;

        % exact derivative of the polynomial above
        theta_rate = w0 + w1 * t_m + w2 * t_m .* t_m;
    end

    traj = struct( ...
        't_m',          t_m, ...
        'theta',        theta, ...
        'theta_rate',   theta_rate, ...
        'w0',           w0, ...
        'w1',           w1, ...
        'w2',           w2, ...
        'T',            T, ...
        'maneuvering',  maneuvering);
end
