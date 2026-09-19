function x_hat = sbl( ...
    y, ...      % the measurement \in \mathbb{C}^{Nmeas}
    A, ...      % the forward operator / sensing matrix (Nmeas x N*N)
    ~, ...      % (was N) grid size; unused -- pixel count is taken from size(A,2)
    varargin ...
    )

    % SBL  Sparse Bayesian Learning (EM / evidence-maximization) for the
    % complex-valued sparse recovery problem  y = A*x + noise.
    %
    % Unlike OMP (greedy, needs sparsity K) and ISTA (convex, needs tau),
    % SBL is a Type-II maximum-likelihood method: it places a zero-mean
    % Gaussian prior x_i ~ CN(0, alpha_i^{-1}) with a SEPARATE precision
    % alpha_i per pixel, then estimates the alpha_i (and the noise
    % variance sigma^2) by maximizing the marginal likelihood (the
    % "evidence"). Sparsity emerges automatically: for irrelevant pixels
    % alpha_i -> infinity, so the posterior mass of x_i collapses onto 0
    % and the column is pruned (Automatic Relevance Determination, ARD).
    % There is thus no sparsity level or regularization weight to set.
    %
    % Seminal / keyed references (equation numbers cited in the steps):
    %   [Tip01] M. E. Tipping, "Sparse Bayesian learning and the relevance
    %          vector machine," J. Mach. Learn. Res., 1:211-244, 2001.
    %          (origin of SBL / the RVM; hierarchical prior + evidence
    %          maximization; the update equations.)
    %   [JXC08] S. Ji, Y. Xue & L. Carin, "Bayesian compressive sensing,"
    %          IEEE Trans. Signal Process., 56(6):2346-2356, 2008.
    %          (SBL specialized to the compressive-sensing model y = A*x;
    %          the notation this file follows most closely.)
    % Supporting references:
    %   [WR04]  D. P. Wipf & B. D. Rao, "Sparse Bayesian learning for basis
    %          selection," IEEE Trans. Signal Process., 52(8):2153-2164,
    %          2004.  (SBL as sparse recovery; why it tolerates dictionary
    %          coherence better than OMP / L1 -- the reason it is worth
    %          adding to the coherence study.)
    %   [Mac92] D. J. C. MacKay, "Bayesian interpolation," Neural Comput.,
    %          4(3):415-447, 1992.  (evidence procedure / the faster
    %          fixed-point hyperparameter update, [Tip01] Eq. (16).)
    %   [TF03]  M. E. Tipping & A. C. Faul, "Fast marginal likelihood
    %          maximisation for sparse Bayesian models," AISTATS, 2003.
    %          (constructive add/delete speed-up; not implemented here.)
    %   [Tib96] R. Tibshirani, "Regression shrinkage and selection via the
    %          lasso," J. R. Stat. Soc. B, 58(1):267-288, 1996.
    %
    % NOTE ON COMPLEX ADAPTATION: [Tip01] and [JXC08] are real-valued.
    % Here x, y, A are complex, so the model is the circularly-symmetric
    % complex Gaussian:  prior x_i ~ CN(0, alpha_i^{-1}), likelihood
    % y ~ CN(A*x, sigma^2 I). All transposes become Hermitian ('), and the
    % real-valued second moment mu_i^2 in the updates becomes |mu_i|^2.
    % The update FORMS are otherwise unchanged (both real and complex
    % circular Gaussians give the EM update alpha_i = 1/E[|x_i|^2]).
    %
    % NOTE ON NOTATION / DIMENSIONS: paper notation clashes with this code
    % (their "N"/"K" is the number of MEASUREMENTS; our N is the image
    % side). Local names used here:
    %   Kpix  = N*N       number of basis functions / pixels (cols of A)
    %                     -> the "N" of [Tip01] Eq. (5), the "N" of [JXC08]
    %   Nmeas = size(A,1) number of measurements (rows of A)
    %                     -> the "N" of [Tip01] Eq. (18), the "K" of [JXC08]
    %
    % NOTE ON COMPUTATION: the posterior covariance Sigma of [Tip01]
    % Eq. (12) / [JXC08] Eq. (10) is inverted directly over the current
    % active set (the |active|-by-|active| Hessian A + beta*As'*As, via
    % Cholesky, an O(M^3) routine, exactly as stated in [Tip01] App. A.1).
    % Pruning columns whose alpha_i has diverged ([Tip01] App. B.1 /
    % Eq. (51); [JXC08] Sec. III-B) shrinks the active set toward the
    % sparsity level, keeping the inverse tractable as iterations proceed.

    if ~isempty(varargin);  max_iter    = varargin{1}; else; max_iter    = 300;  end
    if numel(varargin) > 1; tol         = varargin{2}; else; tol         = 1e-4; end
    if numel(varargin) > 2; prune_ratio = varargin{3}; else; prune_ratio = 1e9;  end

    Kpix  = size(A,2);
    Nmeas = size(A,1);

    % --- Initialization -----------------------------------------------
    % Uniform (a=b=c=d=0) Gamma hyperpriors are assumed, i.e. improper
    % 1/alpha and 1/sigma^2 priors, so no hyperprior constants are set
    % ([Tip01] Sec. 2.3 & Eq. (32); [JXC08] after Eq. (13)). Start every
    % pixel active with unit precision, and set the noise variance to a
    % small fraction of the signal power.
    alpha  = ones(Kpix,1);          % per-pixel precisions alpha_i
    sigma2 = 0.1 * var(y);          % noise variance sigma^2 (= 1/alpha_0)
    active = (1:Kpix).';            % indices of not-yet-pruned columns

    S  = active;                    % active set used for the output
    mu = zeros(0,1);                % posterior mean on the active set

    for iter = 1:max_iter

        S    = active;
        As   = A(:,S);              % Nmeas x |S| active dictionary
        alph = alpha(S);            % |S| x 1 active precisions

        % --- E-step: posterior over the (active) weights --------------
        % Posterior is CN(mu, Sigma) with beta = 1/sigma^2 (= alpha_0):
        %   Sigma = (beta*As'*As + diag(alph))^{-1}   [Tip01 Eq. (12); JXC08 Eq. (10)]
        %   mu    = beta*Sigma*As'*y                  [Tip01 Eq. (13); JXC08 Eq.  (9)]
        beta   = 1 / sigma2;
        Sigma  = inv(beta * (As' * As) + diag(alph));      % |S| x |S|
        mu     = beta * Sigma * (As' * y);                 % posterior mean
        Sig_ii = real(diag(Sigma));                        % diag of Sigma

        % "well-determinedness" gamma_i = 1 - alpha_i*Sigma_ii in [0,1].
        % Both [Tip01] Eq. (16) and [JXC08] Eq. (12) are the fixed-point
        % update alpha_i^new = gamma_i/mu_i^2; gamma_i itself is defined,
        % unnumbered, in the sentence immediately following each equation
        % ("...we have defined the quantities gamma_i by: gamma_i = 1 -
        % alpha_i*Sigma_ii...").
        gamma = 1 - alph .* Sig_ii;

        % --- M-step: re-estimate hyperparameters (EM updates) ---------
        %
        % Precision update. General Gamma(a,b) hyperprior form, obtained
        % by directly differentiating the log-evidence objective (36) and
        % solving for alpha_i:
        %   alpha_i^new = (1+2a) / (mu_i^2 + Sigma_ii + 2b)     [Tip01 Eq. (43)]
        % The text notes this "is equivalent to an expectation-maximisation
        % (EM) update"; the formal EM derivation (treating w as the hidden
        % variable, Appendix A.3) gives the identical formula:
        %   alpha_i = (1+2a) / (<w_i^2> + 2b), <w_i^2> = Sigma_ii+mu_i^2
        %                                                       [Tip01 Eq. (48)]
        % With uniform hyperpriors a=b=0 (Sec. 2.3, "a=b=c=d=0"):
        %   alpha_i^new = 1 / (mu_i^2 + Sigma_ii)
        % (The faster non-EM MacKay fixed point alpha_i=gamma_i/mu_i^2 of
        % [Tip01] Eq. (16) / [JXC08] Eq. (12), computed above as `gamma`,
        % is the alternative -- not used for this update.)
        alpha_new = 1 ./ (abs(mu).^2 + Sig_ii);

        % Noise-variance update. General Gamma(c,d) hyperprior form, via
        % the same formal EM derivation (Appendix A.3, Eq. (49)):
        %   (sigma^2)^new = (||t-Phi*mu||^2 + (sigma^2)^old*sum_i(gamma_i) + 2d)
        %                   / (N + 2c)                          [Tip01 Eq. (50)]
        % With uniform hyperpriors c=d=0:
        %   (sigma^2)^new = (||y-As*mu||^2 + (sigma^2)^old*sum(gamma)) / Nmeas
        % (This is NOT the same formula as the non-EM main-text update
        % (sigma^2)^new = ||t-Phi*mu||^2/(N-sum_i(gamma_i)), [Tip01] Eq. (18)
        % / [JXC08] Eq. (13) -- the EM and non-EM sigma^2 updates differ
        % per-iteration, unlike the alpha case where Eq. (43) and (48)
        % coincide exactly.)
        resid  = y - As * mu;
        sigma2 = (norm(resid)^2 + sigma2 * sum(gamma)) / Nmeas;

        progress_bar('SBL', iter, max_iter);

        % --- Pruning (ARD) --------------------------------------------
        % As alpha_i diverges the pixel is a posteriori certain to be 0,
        % so its column is removed ([Tip01] App. B.1 / Eq. (51); [JXC08]
        % Sec. III-B). We drop columns whose precision has grown many
        % orders of magnitude beyond the most relevant one.
        alpha(S) = alpha_new;
        keep     = alpha_new < prune_ratio * min(alpha_new);
        active   = S(keep);

        if isempty(active)
            break
        end

        % --- Convergence: small change in the log-precisions ----------
        if max(abs(log(alpha_new(keep)) - log(alph(keep)))) < tol
            S  = S(keep);
            mu = mu(keep);
            break
        end
    end

    progress_bar('SBL', 1, 1);
    fprintf('\n');

    % Scatter the active-set posterior mean back into the full image.
    % Pruned pixels stay exactly zero (their alpha_i -> infinity).
    x_hat    = zeros(Kpix,1);
    x_hat(S) = mu;
end
