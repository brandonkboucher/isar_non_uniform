function file = write_scenario_csv(x_hat, alg_names, alg_label, alg_atoms, ...
    alg_build, meta, file)
% WRITE_SCENARIO_CSV  Write one scenario run's summary table as CSV.
%
%   file = WRITE_SCENARIO_CSV(x_hat, alg_names, alg_label, alg_atoms,
%   alg_build, meta) writes the table ISAR_TESTING_V3 prints -- one row per
%   algorithm that ran -- to results/scenario_summary.csv, with the same
%   columns as the sweep's CSV so the two can be read side by side:
%
%     algorithm,run_s,dict_s,total_s,error_m,atoms,matched,missed,false_alarms
%
%   plus the scenario's seed, scatterer count and unambiguous extent, which are
%   the same on every row but make a row self-describing once several runs are
%   concatenated.
%
%   TOTAL_S is RUN_S + DICT_S: the cost of one measurement when the dictionary
%   is built for it. A run time is wall clock and moves between runs, so do not
%   mix a DICT_S from one run with a RUN_S from another.
%
%   META is a struct with fields seed, num_of_scatterers and Wx. FILE overrides
%   the default path.
%
%   See also ISAR_TESTING_V3, WRITE_SWEEP_CSV.

    if nargin < 7 || isempty(file)
        root = fileparts(fileparts(mfilename('fullpath')));
        outdir = fullfile(root, 'results');
        if ~isfolder(outdir)
            mkdir(outdir);
        end
        file = fullfile(outdir, 'scenario_summary.csv');
    end

    fid = fopen(file, 'w');
    if fid < 0
        error('write_scenario_csv:cannotWrite', 'could not open %s', file);
    end
    closer = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, ['algorithm,run_s,dict_s,total_s,error_m,atoms,' ...
        'matched,missed,false_alarms,seed,n_scatterers,Wx_m\n']);

    for ia = 1:numel(alg_names)
        alg = alg_names(ia);
        if ~(isfield(x_hat, alg) && isfield(x_hat.(alg), 'time_s'))
            continue
        end
        res = x_hat.(alg);

        matched = NaN; missed = NaN; false_alarms = NaN;
        if isfield(res, 'pairs'),        matched      = size(res.pairs, 1); end
        if isfield(res, 'missed'),       missed       = numel(res.missed); end
        if isfield(res, 'false_alarms'), false_alarms = numel(res.false_alarms); end

        fprintf(fid, '"%s",%.6g,%.6g,%.6g,%.6g,%d,%d,%d,%d,%d,%d,%.6g\n', ...
            strrep(char(alg_label(ia)), '"', '""'), ...
            res.time_s, alg_build(ia), res.time_s + alg_build(ia), ...
            res.error, alg_atoms(ia), matched, missed, false_alarms, ...
            meta.seed, meta.num_of_scatterers, meta.Wx);
    end
end
