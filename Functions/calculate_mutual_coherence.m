function mu_mat = calculate_mutual_coherence(...
    A,...
    x_array, ...
    y_array, ...
    target_locations)
    
    % normalize the columns
    A_norm = A ./vecnorm(A,2,1);

    % calculate the gram matrix
    G = A_norm' * A_norm;

    % remove the diagonal elements
    G = G - diag(diag(G));

    % find the mutual coherence of A
    mu = max(abs(G(:)));

    % grid columns of A nearest to the true (off-grid) scatterers
    [~, cix] = min(abs(x_array(:) - target_locations(:,1).'), [], 1);
    [~, rix] = min(abs(y_array(:) - target_locations(:,2).'), [], 1);
    idx = (cix(:)-1)*N + rix(:);   % column-major index into the N-by-N grid
    a_vecs = A(:,idx);

    mu_mat = zeros(Ks);
    for i = 1:Ks
        for j = 1:Ks
            if i == j
                mu_mat(i,j) = 1;
                continue
            end
            ai = A(:,i);
            aj = A(:,j);
        
            % evaluate the inner product between atoms
            mu_mat(i,j) = abs(ai'*aj) / (norm(ai)*norm(aj));
        end
    end

end

