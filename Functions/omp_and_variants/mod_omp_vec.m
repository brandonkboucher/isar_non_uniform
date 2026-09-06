
% This is a modified OMP algorithm that identifies and
% removes aliased scatterers from the reconstructed latent
% image. For each selected atom, the algorithm checks the
% correlation of its ambiguous doppelgangers. The
% non-uniform phase progression and the range cell migration
% breaks the ambiguity between the atoms, and the atom which
% correlates best with the residual is selected. If the
% selected atom is identified as ambiguous, its added to the
% sensing matrix and to the latent image.

function x_hat = mod_omp_vec( ...
    y,...       % measurement
    A,...       % gridded sensing matrix
    sparsity, ...      % sparsity
    grid, ...   % grid parameters
    u0, ...     % (m) distance from radar to target rotation axis
    theta_m, ...% [M x 1] yawing angle as function of slow-time
    f_hat_l, ...% [L x 1]
    fc, ...     % (Hz) center frequency
    num_of_amb, ... % the number of ambiguities to check
    options ... % additional scenario options
    )

    % check if the number of ambiguities is odd
    if mod(num_of_amb,2) == 0
        error('The number of ambiguities cannot be even to run mod OMP.')
    end

    % define the dimensions of the matrices
    [ML, K] = size(A);

    % initialize the array containing the approximate non-zero
    % uambiguous indices and track which atoms are
    % unambiguous
    selected_idx = []; 

    % create an array of selected atoms in the unambiguous
    % latent image in order to ensure you aren't
    % re-selecting the same atom twice
    selected_unambiguous_idx = [];

    % define the estimated latent image with corresponding
    % ambiguities i.e. if num_of_amb then x_hat is structure
    % as [-1, 0, 1]
    x_hat = zeros(num_of_amb * K,1);
    
    % create an estimated sensing matrix, composed of
    % selected atoms
    As = zeros(ML, sparsity);
    
    r = y;
    for i = 1:sparsity
    
        % inner product of the measurement and basis vectors
        c = A' * r;
        c(selected_unambiguous_idx) = 0;

        % find the index that maximizes the inner product
        [~, l] = max(abs(c), [] , "all");

        % check the correlation of the selected atom's
        % ambiguous doppelgangers
        [selected_a, selected_ambiguity] ...
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

        if selected_ambiguity ~= 0
            fprintf('Ambiguous scatterer identified.')
        end

        % add new selected atom to sensing matrix
        As(:,i) = selected_a;

        % save index of ambiguous or unambiguous scatterer,
        % assume the num_of_amb is odd, then allows us to
        % shift the index over depending on which ambiguity
        % the scatterer lands in
        selected_idx = [selected_idx; ...
            l + (selected_ambiguity + floor(num_of_amb/2)) * K];
        selected_unambiguous_idx = [selected_unambiguous_idx; l];

        % calculate the residual after adding the new atom
        coef = As(:,1:i) \ y;
        r = y - As(:,1:i) * coef;

        if any(isnan(r))
            error('Mod-OMP has nans.')
        end

        progress_bar('Mod-OMP', i, sparsity);
    end

    % use this basis vector to calculate the least squares
    % solution for the complex reflectivity
    x_hat(selected_idx) = As \ y;

    fprintf('\n');
end

