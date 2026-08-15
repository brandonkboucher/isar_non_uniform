
% This code is based on the paper Newtonized Orthogonal
% Matching Pursuit: Frequency Estimation Over the Continuum.
% The code is structure using a "Detection" step which is
% the standard OMP selection of the atom that maximizes its
% projection onto the residual. The second step "Refinement"
% applies Newton's method to refine the frequency
% estimation.

% Each atom is defined as:
% a = exp(-j * 2 * pi * (fc + \hat{f}) * \tau_{k}(t_{m}))

function [alpha_hat, p_hat] = nomp_vec( ...
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
    fc ...          % center frequency
    )

    % initialize the array containing the approximate non-zero
    % indices
    alpha_hat = zeros(sparsity,1);
    p_hat = zeros(2,sparsity);
    
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

        %---------------- Refinement ----------------

        x0 = xk(k);
        y0 = yk(k);

        % perform Newton's method for optimization
        p_hat(:,iatom) = newton_method(...
            r, x0, y0, u0, theta_m, f_hat_l, fc, Rs);

        % build the dictionary        
        AS = build_AS(...
            p_hat(:,1:iatom), u0, theta_m, f_hat_l, fc);

        % use this basis vector to calculate the least squares
        % solution
        alpha_hat(1:iatom) = AS \ s;

        %-------------- Cyclic Refinement --------------

        for icycle = 1:Rc
            for jatom = 1:iatom
    
                idx = 1:iatom;
                idx(jatom) = [];
    
                % extract position and reflectivity sans the
                % selected atom
                x0 = p_hat(1,jatom); y0 = p_hat(2,jatom);
    
                % build the dictionary        
                AS = build_AS(...
                    p_hat(:,idx), u0, theta_m, f_hat_l, fc);
    
                % form measurement from other atoms
                rj = s - AS * alpha_hat(idx);
    
                % perform Newton's method for optimization
                p_hat(:,jatom) = newton_method(...
                    rj, x0, y0, u0, theta_m, f_hat_l, fc, Rs);
    
            end

            % build the dictionary        
            AS = build_AS(...
                p_hat(:,1:iatom), u0, theta_m, f_hat_l, fc);
    
            % use this basis vector to calculate the least squares
            % solution
            alpha_hat(1:iatom) = AS \ s;

        end
        
        r = s - AS * alpha_hat(1:iatom);

        if any(isnan(r))
            error('OMP has nans.')
        end

        progress_bar('NOMP', iatom, sparsity);
    end
    fprintf('\n');
end

