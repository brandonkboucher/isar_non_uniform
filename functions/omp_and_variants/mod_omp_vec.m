
% This is a modified OMP algorithm that identifies and
% removes aliased scatterers from the reconstructed latent
% image. For each selected atom, the algorithm checks the
% correlation of its ambiguous doppelgangers. The
% non-uniform phase progression and the range cell migration
% breaks the ambiguity between the atoms, and the atom which
% correlates best with the residual is selected. If the
% selected atom is identified as ambiguous, its added to the
% sensing matrix and to the latent image.

function [alpha_hat, p_hat] = mod_omp_vec( ...
    s,...       % measurement
    A,...       % gridded sensing matrix
    sparsity, ...      % sparsity
    grid, ...   % grid parameters
    u0, ...     % (m) distance from radar to target rotation axis
    theta_m, ...% [M x 1] yawing angle as function of slow-time
    f_hat_l, ...% [L x 1]
    fc, ...     % (Hz) center frequency
    num_of_amb, ... % the number of ambiguities to check
    sc, ... % scenario parameters, optimization method
    options ... % additional scenario options
    )

    % check if the number of ambiguities is odd
    if mod(num_of_amb,2) == 0
        error('The number of ambiguities cannot be even to run mod OMP.')
    end
    
    % define the dimensions of the matrices
    [ML, K] = size(A); 

    % create an array of selected atoms in the unambiguous
    % latent image in order to ensure you aren't
    % re-selecting the same atom twice
    selected_unambiguous_idx = [];
    
    % create an estimated sensing matrix, composed of
    % selected atoms
    As = zeros(ML, sparsity);
    
    % define an array containing the approximate positions of 
    % the scatterers, used by PROMP
    p_hat = zeros(2,sparsity);

    % define an array containing the complex reflectivity of
    % each of the scatterers
    alpha_hat = zeros(sparsity,1);

    % Section IV-A stops every algorithm when the signal residual reaches the
    % noise level; Inf when no threshold is set, leaving the sparsity cap
    tau = residual_stop_threshold(options);

    % how many atoms were actually kept, so an early stop returns only those
    n_kept = sparsity;

    r = s;
    for i = 1:sparsity
    
        % inner product of the measurement and basis vectors
        c = A' * r;
        c(selected_unambiguous_idx) = 0;

        % find the index that maximizes the inner product
        [~, l] = max(abs(c), [] , "all");

        % check the correlation of the selected atom's
        % ambiguous doppelgangers
        if strcmpi(sc.optimization_method, 'offset')
            [selected_a, ~, x, y] ...
                = check_ambiguity_of_an_atom(...
                    r, ...          % [M x 1] residual
                    grid.xk(l), ... % (m) crossrange position
                    grid.yk(l), ... % (m) range position
                    grid.Wx, ...    % (m) unambiguous crossrange extent
                    grid.cross_range_pixel_res, ... % (m) crossrange pixel size
                    num_of_amb, ... % the number of ambiguities to check
                    u0, ...         % (m) distance from radar to target rotation axis
                    theta_m, ...    % [M x 1] yawing angle as function of slow-time
                    f_hat_l, ...    % [L x 1]
                    fc, ...         % (Hz) center frequency
                    options ...     % additional scenario options
                    );
        elseif strcmpi(sc.optimization_method, 'newtons')
            [selected_a, ~, x, y] = check_ambiguity_and_apply_newtons(...
                r, ...          % [ML x 1] residual
                grid.xk(l), ... % (m) crossrange position of selected atom
                grid.yk(l), ... % (m) range position of selected atom
                grid.Wx, ...    % (m) unambiguous crossrange extent
                num_of_amb, ... % the number of ambiguities to check
                u0, ...         % (m) distance from radar to target rotation axis
                theta_m, ...    % [M x 1] yawing angle as function of slow-time
                f_hat_l, ...    % [L x 1]
                fc, ...         % (Hz) center frequency
                sc.num_optimization_steps, ... % number of refinement steps
                options ...     % additional scenario options
                );
        elseif strcmpi(sc.optimization_method, 'newtons_and_offset')
            [selected_a, ~, x, y] = check_ambiguity_and_apply_newtons_with_offsets(...
                r, ...          % [ML x 1] residual
                grid.xk(l), ... % (m) crossrange position of selected atom
                grid.yk(l), ... % (m) range position of selected atom
                grid.Wx, ...    % (m) unambiguous crossrange extent
                num_of_amb, ... % the number of ambiguities to check
                u0, ...         % (m) distance from radar to target rotation axis
                theta_m, ...    % [M x 1] yawing angle as function of slow-time
                f_hat_l, ...    % [L x 1]
                fc, ...         % (Hz) center frequency
                sc.num_optimization_steps, ... % number of refinement steps
                grid.cross_range_pixel_res, ... % (m) crossrange pixel size
                sc.num_offsets_pixels, ... % number of pixels to traverse using Newton's method
                options ...     % additional scenario options
                );
        else
            error('You must select either newtons, offset or newtons_and_offset to execute mod_omp')
        end
        p_hat(:,i) = [x;y];

        % add new selected atom to sensing matrix
        As(:,i) = selected_a;

        % save index of ambiguous or unambiguous scatterer,
        % assume the num_of_amb is odd, then allows us to
        % shift the index over depending on which ambiguity
        % the scatterer lands in
        selected_unambiguous_idx ...
            = [selected_unambiguous_idx; l];

        % calculate the residual after adding the new atom
        alpha_hat = As(:,1:i) \ s;
        r = s - As(:,1:i) * alpha_hat;

        if any(isnan(r))
            error('Mod-OMP has nans.')
        end

        if options.debug_printing
            progress_bar('Mod-OMP', i, sparsity);
        end

        % the residual has reached the noise level
        if norm(r) <= tau
            n_kept = i;
            break
        end
    end

    % an early stop leaves the tail of the preallocated arrays untouched
    p_hat     = p_hat(:, 1:n_kept);
    alpha_hat = alpha_hat(1:n_kept);

    if options.debug_printing
        fprintf('\n');
    end
end

