function x_hat = omp_vec( ...
    y, ... % the measurement \in \mathbb{C}^{M}
    A, ... % the forward operator function
    K ... % Sparsity
    )

    % initialize the array containing the approximate non-zero
    % indices
    Lambda = []; 
    x_hat = zeros(size(A,2),1);
    
    % for now, assume we have knowledge of the number of
    % nonzeros
    r = y;
    for i = 1:K
    
        % inner product of the measurement and basis vectors
        c = A' * r;
    
        % find the index that maximizes the inner product
        [~, l] = max(abs(c), [] , "all");
        Lambda = [Lambda; l];
    
        % 3) Build restricted sensing matrix A_S (M x |Lambda|)        
        AS = zeros(size(y,1), length(Lambda));        
        for j = 1:length(Lambda)            
            ej = zeros(size(A,2),1);            
            ej(Lambda(j)) = 1;            
            AS(:,j) = A * ej;       % forward-proj an impulse        
        end

        % use this basis vector to calculate the least squares
        % solution
        x_hat(Lambda) = inv(AS'*AS)*AS'*y;
        r = y - A * x_hat;

        if any(isnan(r))
            error('OMP has nans.')
        end

        progress_bar('OMP', i, K);
    end
    fprintf('\n');
end

