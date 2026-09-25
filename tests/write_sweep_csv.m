function files = write_sweep_csv(res, outdir)
% WRITE_SWEEP_CSV  Write the paper sweep's results as CSV.
%
%   files = WRITE_SWEEP_CSV(res, outdir) takes the struct RUN_PAPER_SWEEP
%   returns and writes two files into OUTDIR (default: results/):
%
%     paper_sweep_summary.csv   one row per arm -- the table the paper quotes:
%                               median, mean over the draws below 1 m, 90th
%                               percentile, gross-failure count, run time per
%                               draw, dictionary build time, per-measurement
%                               cost (dictionary + one run) and dictionary size
%
%     paper_sweep_errors.csv    one row per draw: the draw index, its seed, and
%                               every arm's RMS position error, so the
%                               distributions can be replotted without MATLAB
%
%   A run of other than 500 draws (the paper's sweep size) is written to
%   paper_sweep_summary_<n>draws.csv instead, so a smoke test cannot overwrite
%   the full run's files.
%
%   Written with fprintf rather than writetable so no toolbox is needed, and
%   with %.6g so a value read back matches the .mat baseline to the tolerance
%   TEST_PAPER_SWEEP asserts.
%
%   See also RUN_PAPER_SWEEP, TEST_PAPER_SWEEP.

    if nargin < 2 || isempty(outdir)
        root   = fileparts(fileparts(mfilename('fullpath')));
        outdir = fullfile(root, 'results');
    end
    if ~isfolder(outdir)
        mkdir(outdir);
    end

    % a partial run gets its own filenames, so a 20-draw smoke test cannot
    % overwrite the full sweep's files
    tag = '';
    if res.n_iter ~= 500
        tag = sprintf('_%ddraws', res.n_iter);
    end

    gross = 1.0;   % as in test_paper_sweep: a scatterer in the wrong place
    td     = [res.time_A_mod, res.time_A_base];
    run_s  = res.T / res.n_iter;
    dict_s = td(res.arm_dict);
    atoms  = [res.K_mod, res.K_base];
    atoms  = atoms(res.arm_dict);

    % ---- per-arm summary ----------------------------------------------
    summary_file = fullfile(outdir, ['paper_sweep_summary' tag '.csv']);
    fid = fopen(summary_file, 'w');
    if fid < 0
        error('write_sweep_csv:cannotWrite', 'could not open %s', summary_file);
    end
    closer = onCleanup(@() fclose(fid));
    fprintf(fid, ['arm,median_m,mean_below_1m_m,p90_m,gross_draws,' ...
        'run_ms_per_draw,dict_s,per_measurement_s,atoms,n_draws\n']);
    for a = 1:numel(res.labels)
        e = res.E(:,a);
        fprintf(fid, '%s,%.6g,%.6g,%.6g,%d,%.6g,%.6g,%.6g,%d,%d\n', ...
            csv_field(res.labels(a)), median(e), mean(e(e < gross)), pct90(e), ...
            sum(e >= gross), 1e3*run_s(a), dict_s(a), run_s(a) + dict_s(a), ...
            atoms(a), res.n_iter);
    end
    clear closer

    % ---- per-draw errors ----------------------------------------------
    errors_file = fullfile(outdir, ['paper_sweep_errors' tag '.csv']);
    fid = fopen(errors_file, 'w');
    if fid < 0
        error('write_sweep_csv:cannotWrite', 'could not open %s', errors_file);
    end
    closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, 'draw,seed');
    for a = 1:numel(res.labels)
        fprintf(fid, ',%s', csv_field(res.labels(a)));
    end
    fprintf(fid, '\n');
    for it = 1:res.n_iter
        fprintf(fid, '%d,%d', it, res.base_seed + it - 1);
        fprintf(fid, ',%.6g', res.E(it,:));
        fprintf(fid, '\n');
    end

    files = {summary_file, errors_file};
end

% ------------------------------------------------------------------------
function s = csv_field(label)
% quote a field, and double any quote inside it, so arm names holding commas
% ("mod-OMP orig, new offsets") survive a round trip
    s = char(label);
    s = ['"' strrep(s, '"', '""') '"'];
end

% ------------------------------------------------------------------------
function p = pct90(e)
% 90th percentile by nearest rank; avoids prctile and the Statistics toolbox
    e = sort(e(~isnan(e)));
    if isempty(e), p = NaN; return, end
    p = e(max(1, ceil(0.90 * numel(e))));
end
