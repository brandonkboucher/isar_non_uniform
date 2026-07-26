function x_hat = lasso_ista( ...
    y, ...      % the measurement \in \mathbb{C}^{M}
    A, ...      % the forward operator / sensing matrix
    ~, ...      % (was N) grid size; unused -- pixel count is taken from size(A,2)
    tau, ...    % L1 regularization weight (sparsity-promoting)
    varargin ...
    )

    % LASSO_ISTA  Iterative Soft-Thresholding Algorithm for the
    % complex-valued LASSO / basis-pursuit denoising (BPDN) problem
    %
    %   x_hat = argmin_x  0.5*||y - A*x||_2^2 + tau*||x||_1        (objective)
    %
    % This is the constant-stepsize proximal-gradient (forward-backward
    % splitting) method: a gradient (forward) step on the smooth
    % data-fidelity term f(x) = 0.5*||y - A*x||_2^2 followed by a
    % proximal (backward) step on g(x) = tau*||x||_1, whose proximal
    % operator is complex soft-thresholding / shrinkage.
    %
    % The step comments below cite equations of the two seminal papers:
    %   [BT09] A. Beck & M. Teboulle, "A fast iterative shrinkage-
    %          thresholding algorithm for linear inverse problems,"
    %          SIAM J. Imaging Sci., 2(1):183-202, 2009.
    %   [DDD04] I. Daubechies, M. Defrise & C. De Mol, "An iterative
    %          thresholding algorithm for linear inverse problems with a
    %          sparsity constraint," Comm. Pure Appl. Math.,
    %          57(11):1413-1457, 2004.  (the seminal ISTA convergence
    %          result; equation numbers cited here are from the arXiv
    %          version math/0307152v2.)
    % Supporting references:
    %   [CW05]  P. L. Combettes & V. R. Wajs, "Signal recovery by proximal
    %          forward-backward splitting," Multiscale Model. Simul.,
    %          4(4):1168-1200, 2005.  (proximal forward-backward framework)
    %   [Tib96] R. Tibshirani, "Regression shrinkage and selection via the
    %          lasso," J. R. Stat. Soc. B, 58(1):267-288, 1996.  (LASSO
    %          objective; = BPDN of [CDS98])
    %   [CDS98] S. S. Chen, D. L. Donoho & M. A. Saunders, "Atomic
    %          decomposition by basis pursuit," SIAM J. Sci. Comput.,
    %          20(1):33-61, 1998.  (basis pursuit / BPDN)
    %   [Don95] D. L. Donoho, "De-noising by soft-thresholding," IEEE
    %          Trans. Inf. Theory, 41(3):613-627, 1995.  (soft-threshold)
    %
    % NOTE ON NORMALIZATION: this code uses the 1/2-scaled data term
    % f(x) = 0.5*||y - A*x||_2^2, so grad f = A'*(A*x - y) and the
    % Lipschitz constant is L = lambda_max(A'*A). [BT09] instead uses
    % f(x) = ||A*x - b||^2 (no 1/2), so its gradient carries a factor 2,
    % i.e. grad f = 2*A'*(A*x - b) in Eq. (1.4), and its Lipschitz
    % constant is L(f) = 2*lambda_max(A'*A) (Example 2.2). The two
    % formulations produce identical iterates after rescaling; the code
    % below is the 1/2-normalized form (the convention of [CW05]).
    %
    % The [BT09] objective is Eq. (1.3): min ||A*x - b||^2 + lambda*||x||_1,
    % with their lambda playing the role of tau here.  The [DDD04] form is
    % Eq. (8): ||K*z - v||^2 + mu*||z||_{l1}.

    if ~isempty(varargin); max_iter = varargin{1}; else; max_iter = 200;  end
    if numel(varargin) > 1; tol      = varargin{2}; else; tol      = 1e-6; end

    K = size(A,2);

    % --- Stepsize selection -------------------------------------------
    % ISTA with constant stepsize uses t = 1/L, where L is a Lipschitz
    % constant of grad f. For f(x) = 0.5*||y - A*x||^2 this is
    % L = lambda_max(A'*A) (cf. [BT09] Example 2.2, which gives
    % L(f) = 2*lambda_max(A'*A) for their unscaled f). The convergence
    % condition on the stepsize is t in (0, 1/||A'*A||), stated in [BT09]
    % Sec. 2.1 (just after Example 2.2). We obtain L = lambda_max(A'*A)
    % by power iteration on the operator A'*A so the K-by-K matrix A'*A
    % is never formed.
    v = randn(K,1) + 1j*randn(K,1);
    v = v / norm(v);
    for i = 1:20
        v = A' * (A * v);
        Lip = norm(v);              % Rayleigh-quotient estimate of lambda_max
        v = v / Lip;
    end
    step = 1 / Lip;                 % t = 1/L

    % --- Step 0 [BT09, "ISTA with constant stepsize" box]: init x_0 ----
    x_hat = zeros(K,1);

    for iter = 1:max_iter

        % --- Step k [BT09 Eq. (1.4) / Eq. (3.1); DDD04 Eq. (19)] -------
        % ISTA update: x_k = T_{tau*t}( x_{k-1} - t*grad f(x_{k-1}) ).
        % [BT09] Eq. (3.1) is x_k = p_L(x_{k-1}); the text after it notes
        % that for f=||Ax-b||^2, g=lambda||x||_1 this reduces to the
        % shrinkage form Eq. (1.4) with t = 1/L(f). [DDD04] Corollary 2.3
        % Eq. (19): f^n = S_{w,p}(f^{n-1} + K'(g - K f^{n-1})).
        %
        % (a) gradient (forward) step on f.  grad f(x) = A'*(A*x - y);
        %     this is the A'*(A*x - b) factor inside [BT09] Eq. (1.4)
        %     (their 2*A'*(A*x-b) under the unscaled-f convention), and
        %     equals the -K'(K f - g) direction of [DDD04] Eq. (19).
        grad  = A' * (A * x_hat - y);
        % (b) proximal (backward) step: soft-threshold at level tau*t.
        %     The prox of tau*||.||_1 is soft-thresholding, [BT09] Eq.
        %     (1.5) (real case) / [DDD04] Eqs. (5),(16); see soft_threshold.
        x_new = soft_threshold(x_hat - step * grad, tau * step);

        progress_bar('LASSO', iter, max_iter);

        % relative-change stopping criterion (function values are
        % monotone non-increasing for ISTA, [BT09] Remark 3.1).
        if norm(x_new - x_hat) / max(norm(x_hat), eps) < tol
            x_hat = x_new;
            break
        end
        x_hat = x_new;
    end
    progress_bar('LASSO', 1, 1);
    fprintf('\n');
end

function x_out = soft_threshold(v, thresh)
    % Complex soft-thresholding (shrinkage) operator, the proximal
    % operator of thresh*||.||_1.
    %
    % Real scalar form, T_alpha(x) = (|x| - alpha)_+ sgn(x):
    %   [BT09] Eq. (1.5); [DDD04] Eqs. (5) and (16); [Don95].
    % Complex form used here, shrink magnitude and preserve phase:
    %   T_thresh(v) = max(|v| - thresh, 0) .* v./|v|
    %   = [DDD04] Remark 2.5:  S_mu(r*e^{i*theta}) = S_mu(r)*e^{i*theta}.
    % (This is the exact prox of the l1 norm on C^K; the real Eq. (1.5)
    % of [BT09] is the special case of real-valued v.)
    mag   = abs(v);
    x_out = max(mag - thresh, 0) .* (v ./ max(mag, eps));
end
