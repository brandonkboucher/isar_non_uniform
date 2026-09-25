function report_sweep(res, gross)
% REPORT_SWEEP  Print the paper sweep's summary table.
%
%   REPORT_SWEEP(res, gross) prints one row per arm of the struct
%   RUN_PAPER_SWEEP returns: median, mean over the draws below GROSS, 90th
%   percentile, gross-failure count, run time per draw and the per-measurement
%   cost (dictionary build + one run). GROSS defaults to 1 m, the point at
%   which a scatterer is in the wrong place rather than imprecisely refined.
%
%   See also RUN_PAPER_SWEEP, WRITE_SWEEP_CSV, TEST_PAPER_SWEEP.

    if nargin < 2 || isempty(gross), gross = 1.0; end

    td  = [res.time_A_mod, res.time_A_base];
    run_s  = res.T / res.n_iter;
    dict_s = td(res.arm_dict);

    fprintf('\n  Wx %.3f m, 2D range offset %.3f m, miss penalty %.3f m\n', ...
        res.Wx, res.range_offset, res.miss_penalty);
    fprintf('  dictionaries: A_mod %d atoms in %.2f s, A_base %d atoms in %.2f s\n\n', ...
        res.K_mod, res.time_A_mod, res.K_base, res.time_A_base);

    fprintf('  %-26s %8s %8s %8s %6s %9s %11s\n', 'arm', 'median', 'mean<1m', ...
        'p90', 'gross', 'run/draw', 'per meas.');
    fprintf('  %s\n', repmat('-', 1, 84));
    for a = 1:numel(res.labels)
        e = res.E(:,a);
        fprintf('  %-26s %8.4f %8.4f %8.4f %6d %7.1f ms %9.3f s\n', res.labels(a), ...
            median(e), mean(e(e < gross)), pct90(e), sum(e >= gross), ...
            1e3*run_s(a), run_s(a) + dict_s(a));
    end
    fprintf('  %s\n', repmat('-', 1, 84));
    fprintf(['  gross = draws with a scatterer in the wrong place (>= %.1f m RMS).\n' ...
        '  per meas. = dictionary build + one run, the cost when the dictionary\n' ...
        '  cannot be reused across measurements.\n'], gross);
end

% ------------------------------------------------------------------------
function p = pct90(e)
% 90th percentile by nearest rank, so the report does not need prctile (and
% with it the Statistics toolbox) -- the same reason
% calculate_reconstruction_error builds its cost matrix by hand
    e = sort(e(~isnan(e)));
    if isempty(e), p = NaN; return, end
    p = e(max(1, ceil(0.90 * numel(e))));
end
