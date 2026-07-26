function x_hat = reconstruct_latent_image(...
    y, a, N, Ks,options)

    % N is the latent-image grid size: scalar for a square N-by-N grid, or
    % [Ny Nx] when the range and crossrange pixel counts differ (which they do
    % once the crossrange extent is set by the ambiguity count).
    if isscalar(N)
        Ny = N; Nx = N;
    else
        Ny = N(1); Nx = N(2);
    end

    if options.execute_omp
        x_hat_omp = omp_vec(y,a,N,Ks);
        x_hat.omp = reshape(x_hat_omp,Ny,Nx);
    end

    % LASSO reconstruction via basic ISTA
    if options.execute_itsa
        lasso_tau = 0.1 * max(abs(a'*y));
        x_hat_lasso = lasso_ista(y,a,N,lasso_tau);
        x_hat.lasso = reshape(x_hat_lasso,Ny,Nx);
    end

    % Sparse Bayesian Learning (EM-SBL) reconstruction
    if options.execute_sbl
        x_hat_sbl = sbl(y,a,N);
        x_hat.sbl = reshape(x_hat_sbl,Ny,Nx);
    end

    % filtered backprojection using pseudo-inverse
    if options.execute_fbp
        x_hat_bp = a'*y;
        x_hat.fbp = reshape(x_hat_bp,Ny,Nx);
    end

end

