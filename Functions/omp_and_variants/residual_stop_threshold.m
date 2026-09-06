function tau = residual_stop_threshold(options)
% RESIDUAL_STOP_THRESHOLD  The ||r||_2 at which a greedy pursuit stops.
%
%   tau = RESIDUAL_STOP_THRESHOLD(options) returns options.residual_threshold
%   when it is set, and -Inf otherwise. The pursuits test  norm(r) <= tau , so
%   -Inf is the value that never fires and leaves the sparsity cap governing.
%   (Inf would be exactly wrong: it stops after the first atom.)
%
%   Nguyen et al. (IEEE T-AES 55(6), 2019) stop every algorithm in their
%   comparison "when the signal residual reaches the noise level" (Section
%   IV-A) rather than at a fixed atom count. That is the criterion this
%   threshold carries: for complex noise of per-sample variance sigma^2 over
%   N samples, E{||e||^2} = N*sigma^2, so the caller sets
%
%       options.residual_threshold = sqrt(N * sigma2).
%
%   The criterion is only meaningful when the measurement actually has noise.
%   With noiseless data there is no floor for the residual to reach, so the
%   caller leaves the field empty and the sparsity cap governs.
%
%   See also OMP_VEC, PROMP_VEC, NOMP_VEC, MOD_OMP_VEC.

    tau = -inf;

    if nargin >= 1 && isstruct(options) ...
            && isfield(options, 'residual_threshold') ...
            && ~isempty(options.residual_threshold) ...
            && isfinite(options.residual_threshold)

        tau = options.residual_threshold;
    end
end
