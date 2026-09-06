
% This code is based on the paper Newtonized Orthogonal
% Matching Pursuit: Frequency Estimation Over the Continuum.
% The code is structure using a "Detection" step which is
% the standard OMP selection of the atom that maximizes its
% projection onto the residual. The second step "Refinement"
% applies Newton's method to refine the frequency
% estimation.

% Each atom is defined as:
% a = exp(-j * 2 * pi * (fc + \hat{f}) * \tau_{k}(t_{m}))

function [alpha_hat, p_hat, p_hat_hist] = nomp_vec( ...
    s, ...          % measurement [ML x 1]
    A, ...          % sensing matrix [ML x K]
    Rs, ...         % number of Newton steps per atom selection
    Rc, ...         % number of cyclic refinements
    sparsity, ...   % sparsity
    xk, ...         % x position for each k-index [K]
    yk, ...         % y position for each k-index [K]
    u0, ...         % center of rotation
    theta_m, ...    % yaw angle as a function of time [M]
    f_hat_l, ...    % range-frequencies [L]
    fc, ...          % center frequency
    options ...
    )

    % initialize the array containing the approximate non-zero
    % indices
    alpha_hat = zeros(sparsity,1);
    p_hat = zeros(2,sparsity);

    save_histories = false; p_hat_hist = [];
    if isfield(options, 'save_histories') ...
            && options.save_histories
        
        save_histories = true;

        atom_save_idx = 1;
        if isfield(options, 'save_atom_idx')
            atom_save_idx = options.save_atom_idx;
        end

        % calculate the total number of steps to save
        number_of_total_steps = Rs + 1 + Rs*Rc * atom_save_idx;

        % keep track of every iteration's position
        p_hat_hist = zeros(2, number_of_total_steps);
    end
      
    % for now, assume we have knowledge of the number of
    % nonzeros
    r = s;
    istep = 1;
    for iatom = 1:sparsity
    
        %---------------- Detection ----------------
        % inner product of the measurement and atoms
        c = A' * r;
    
        % find the index that maximizes its projection onto
        % the residual
        [~, k] = max(abs(c), [] , "all");

        %---------------- Refinement ----------------

        x0 = xk(k);
        y0 = yk(k);

        % perform Newton's method for optimization
        [p_hat(:,iatom), p_hat_ns] = newton_method(...
            r, x0, y0,  u0, theta_m, f_hat_l, fc, Rs, options.use_range_approx);

        if save_histories && iatom == atom_save_idx

            % the grid node the detection step picked, then however many
            % Newton steps were actually accepted (newton_method returns a
            % variable number, so the width is read off the result)
            p_hat_hist(:, istep) = [x0; y0];
            istep = istep + 1;

            ns = size(p_hat_ns, 2);
            p_hat_hist(:, istep:istep+ns-1) = p_hat_ns;
            istep = istep + ns;
        end

        % build the dictionary        
        AS = build_AS(...
            p_hat(:,1:iatom), u0, theta_m, f_hat_l, fc, options.use_range_approx);

        % use this basis vector to calculate the least squares
        % solution
        alpha_hat(1:iatom) = AS \ s;

        %-------------- Cyclic Refinement --------------

        if sparsity ~= 1
            for icycle = 1:Rc
                for jatom = 1:iatom
        
                    idx = 1:iatom;
                    idx(jatom) = [];
        
                    % extract position and reflectivity sans the
                    % selected atom
                    x0 = p_hat(1,jatom); y0 = p_hat(2,jatom);
        
                    % build the dictionary        
                    AS = build_AS(...
                        p_hat(:,idx), u0, theta_m, f_hat_l, fc, options.use_range_approx);
        
                    % form measurement from other atoms
                    rj = s - AS * alpha_hat(idx);
        
                    % perform Newton's method for optimization
                    [p_hat(:,jatom), p_hat_ns] = newton_method(...
                        rj, x0, y0, u0, theta_m, f_hat_l, fc, Rs, options.use_range_approx);
        
                    if save_histories && jatom == atom_save_idx
                        ns = size(p_hat_ns, 2);
                        p_hat_hist(:, istep:istep+ns-1) = p_hat_ns;
                        istep = istep + ns;
                    end

                end
    
                % build the dictionary        
                AS = build_AS(...
                    p_hat(:,1:iatom), u0, theta_m, f_hat_l, fc, options.use_range_approx);
        
                % use this basis vector to calculate the least squares
                % solution
                alpha_hat(1:iatom) = AS \ s;
    
            end
        end
        r = s - AS * alpha_hat(1:iatom);

        if any(isnan(r))
            error('OMP has nans.')
        end

        progress_bar('NOMP', iatom, sparsity);
    end
    fprintf('\n');

    % number_of_total_steps is only an estimate -- rejected Newton steps make
    % the history shorter, more cyclic passes make it longer -- so trim to
    % what was written
    if save_histories
        p_hat_hist = p_hat_hist(:, 1:istep-1);
    end
end

