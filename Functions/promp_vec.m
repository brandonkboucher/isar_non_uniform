
% This code is based on the paper Parameter-Refined OMP for Compressive
% Radar Imaging of Rotating Targets (Nguyen, Dogancay, Tran and Berry,
% IEEE T-AES 55(6), 2019), Table II. The structure is a "Detection" step,
% which is the standard OMP selection of the atom that maximizes its
% projection onto the residual, followed by a joint non-linear least
% squares refinement of the positions AND reflectivities of every atom in
% the current support -- which is what replaces NOMP's per-atom Newton
% step and its cyclic refinement.

% Each atom is defined as:
% a = exp(-j * 2 * pi * (fc + \hat{f}) * \tau_{k}(t_{m}))

function [alpha_hat, p_hat] = promp_vec( ...
    s, ...          % measurement [ML x 1]
    A, ...          % sensing matrix [ML x K]
    Rs, ...         % number of Gauss-Newton steps per atom selection
    Rc, ...         % unused by PROMP (the joint NLS replaces cyclic refinement)
    sparsity, ...   % sparsity
    xk, ...         % x position for each k-index [K]
    yk, ...         % y position for each k-index [K]
    u0, ...         % center of rotation
    theta_m, ...    % yaw angle as a function of time [M]
    f_hat_l, ...    % range-frequencies [L]
    fc ...          % center frequency
    )

    const = Constants;

    % the constraint of eq 20 is set to the conventional radar resolution. the
    % paper uses lambda/2 because for its single-frequency CW geometry that IS
    % the resolution; here range resolution is set by bandwidth and crossrange
    % resolution by the angular span, and the two differ by several times, so
    % the bound is applied per axis. a bound tighter than the resolution
    % truncates any scatterer that sits further than the bound from its grid
    % node, which shows up as an error floor that no amount of SNR removes.
    range_resolution = const.c / (2 * (max(f_hat_l) - min(f_hat_l)));
    sin_span = max(sin(theta_m)) - min(sin(theta_m));
    cross_range_resolution = const.c / (2 * (fc + max(f_hat_l)) * sin_span);
    zeta = [cross_range_resolution; range_resolution];

    % initialize the array containing the approximate non-zero
    % indices
    alpha_hat = zeros(sparsity,1);
    p_hat = zeros(2,sparsity);

    % the constraint of eq 20 is measured from the grid node that was
    % selected for each atom, so those nodes are recorded separately and
    % never updated by the refinement
    p_anchor = zeros(2,sparsity);

    % support set, used to stop the iteration if an atom is reselected
    Lambda = [];
    
    % for now, assume we have knowledge of the number of
    % nonzeros
    r = s;
    for iatom = 1:sparsity
    
        %---------------- Detection ----------------
        % inner product of the measurement and atoms
        c = A' * r;
    
        % find the index that maximizes its projection onto
        % the residual
        [~, k] = max(abs(c), [] , "all");

        % Table II quits the iteration if the selected atom is already in
        % the support, which a dense oversampled grid makes possible
        if ismember(k, Lambda)
            break
        end
        Lambda = [Lambda; k];

        %---------------- Refinement ----------------

        % extract the coarse position estimate
        x0 = xk(k);
        y0 = yk(k);
        p_anchor(:,iatom) = [x0; y0];

        % seed the new atom's reflectivity from the standard least squares
        % projection onto the residual (eq 24a). the paper's atoms are
        % energy-normalised so that this is simply c(k); ours are not, so
        % divide through by the atom energy
        a = compute_atom(x0, y0, u0, theta_m, f_hat_l, fc, false);
        alpha_hat(iatom) = c(k) / (a' * a);

        % the refinement uses a non-linear least squares solver NLS, which
        % jointly refines every atom in the support against the measurement
        [p_hat, alpha_hat] = promp_nls_solver(...
            s, ...          % measurement
            p_hat, ...      % previous estimates (i-1-th)
            x0, ...         % new coarse estimate in x-direction
            y0, ...         % new coarse estimate in y-direction
            alpha_hat, ...  % previous reflectivity estimate
            iatom, ...      % estimate index 
            Rs, ...         % max number of Gauss-Newton iterations
            u0, ...         % [m] distance from radar to scene center
            theta_m, ...    % [rad] angle as a function of slow-time
            f_hat_l, ...    % [Hz] range-frequency
            fc, ...         % [Hz] radar center frequency
            zeta, ...                   % constraint, the conventional resolution per axis
            p_anchor(:,1:iatom) ...     % grid nodes the constraint is measured from
            );

        % build the dictionary        
        AS = build_AS(...
            p_hat(:,1:iatom), u0, theta_m, f_hat_l, fc);

        r = s - AS * alpha_hat(1:iatom);

        if any(isnan(r))
            error('OMP has nans.')
        end

        progress_bar('PROMP', iatom, sparsity);
    end
    fprintf('\n');
end
