function x_hat = omp( ...
    y, ... % the measurement \in \mathbb{C}^{M}
    A, ... % the forward operator function
    AH, ... % the adjoint operator function
    Phi, ... % sampling matrix \Phi \in {0,1}^{M \times N^2}
    N, ... % Dimensionality of latent image (x \in \mathbb{C}^{N \times N})
    K ... % Sparsity
    )

    % initialize the array containing the approximate non-zero
    % indices
    Lambda = []; 
    x_hat = zeros(N,N);
    
    % for now, assume we have knowledge of the number of
    % nonzeros
    r = y;
    for i = 1:K
    
        % inner product of the measurement and basis vectors
        c = AH(r,Phi,N);
    
        % find the index that maximizes the inner product
        [~, l] = max(abs(c), [] , "all");
        Lambda = [Lambda; l];
    
        % 3) Build restricted sensing matrix A_S (M x |Lambda|)        
        AS = zeros(size(y,1), length(Lambda));        
        for j = 1:length(Lambda)            
            ej = zeros(N,N);            
            ej(Lambda(j)) = 1;            
            AS(:,j) = A(ej, Phi);       % forward-proj an impulse        
        end

        % use this basis vector to calculate the least squares
        % solution
        x_hat(Lambda) = inv(AS'*AS)*AS'*y;
        r = y - A(x_hat, Phi);
    end
end

