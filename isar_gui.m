function isar_gui()
% ISAR_GUI  Front end for the ISAR non-uniform sampling reconstruction study.
%
%   ISAR_GUI opens a window whose dropdowns mirror the scenario columns of
%   documentation/scenarios_new.xlsx. Pick a configuration -- or load one of
%   the numbered cases from the spreadsheet -- press Run, and the formed
%   image and the per-scatterer position errors appear for every image
%   formation algorithm you selected.
%
%   The simulation itself is isar_run_scenario, the same function the batch
%   script isar_testing_v2.m drives, so the GUI and the batch sweep produce
%   identical numbers for identical settings.
%
%   Every run is appended to the History tab, so stepping a single dropdown
%   -- Critical to Oversampled, say -- builds a side-by-side comparison you
%   can export to CSV.

    root = fileparts(mfilename('fullpath'));
    addpath(fullfile(root, 'functions'));

    % ---- shared state (visible to all nested functions) ----------------
    S.root      = root;
    S.scenarios = load_scenarios(root);
    S.last      = [];      % output of the most recent run
    S.history   = {};      % one row per run, for the History tab
    S.axesList  = gobjects(0);

    build_ui();
    apply_preset_items();
    sync_enable();
    render_trajectory();

%% ------------------------------------------------------------------ UI

    function build_ui()

        S.fig = uifigure('Name', 'ISAR Non-Uniform Sampling Simulator', ...
            'Position', [80 80 1440 900]);

        outer = uigridlayout(S.fig, [1 2]);
        outer.ColumnWidth = {360, '1x'};
        outer.RowHeight   = {'1x'};

        % ---------------- left: controls --------------------------------
        % every panel sizes to its contents and the column scrolls, so the
        % solver panel is never clipped on a short screen
        left = uigridlayout(uipanel(outer, 'BorderType', 'none'), [4 1]);
        left.RowHeight   = {'fit', 'fit', 'fit', 'fit'};
        left.ColumnWidth = {'1x'};
        left.Scrollable  = 'on';

        build_scenario_panel(left);
        build_algorithm_panel(left);
        build_advanced_panel(left);
        build_run_controls(left);

        % ---------------- right: results --------------------------------
        S.tabs = uitabgroup(outer);

        S.imageTab = uitab(S.tabs, 'Title', 'Formed images');
        S.imageGrid = uigridlayout(S.imageTab, [1 1]);
        placeholder(S.imageGrid, 'Press Run to form an image.');

        S.trajTab = uitab(S.tabs, 'Title', 'Trajectory');
        build_trajectory_tab();

        S.errorTab = uitab(S.tabs, 'Title', 'Position errors');
        build_error_tab();

        S.historyTab = uitab(S.tabs, 'Title', 'History');
        build_history_tab();

        S.logTab = uitab(S.tabs, 'Title', 'Log');
        logGrid = uigridlayout(S.logTab, [1 1]);
        S.logArea = uitextarea(logGrid, 'Editable', 'off', ...
            'FontName', 'Menlo', 'Value', {''});
    end

    function build_scenario_panel(parent)

        p = uipanel(parent, 'Title', 'Scenario');
        g = uigridlayout(p, [9 2]);
        g.ColumnWidth = {170, '1x'};
        g.RowHeight   = repmat({26}, 1, 9);

        % preset picker
        uilabel(g, 'Text', 'Spreadsheet case');
        S.presetDrop = uidropdown(g, 'ValueChangedFcn', @on_preset);

        % the scenario columns themselves
        uilabel(g, 'Text', '# ambiguities w/ scatterers');
        S.ambScatDrop = uidropdown(g, 'Items', {'1', '3'}, ...
            'ItemsData', [1 3], 'ValueChangedFcn', @on_scenario_change);

        uilabel(g, 'Text', '# ambiguities in image former');
        S.ambIfDrop = uidropdown(g, 'Items', {'1', '3'}, ...
            'ItemsData', [1 3], 'ValueChangedFcn', @on_scenario_change);

        uilabel(g, 'Text', 'Scatterer locations');
        S.scatLocDrop = uidropdown(g, 'Items', {'On-Grid', 'Off-Grid'}, ...
            'ValueChangedFcn', @on_scenario_change);

        uilabel(g, 'Text', 'Image former grid density');
        S.gridDensityDrop = uidropdown(g, 'Items', {'Critical', 'Oversampled'}, ...
            'ValueChangedFcn', @on_scenario_change);

        uilabel(g, 'Text', 'Angle rate');
        S.angleRateDrop = uidropdown(g, 'Items', {'Constant', 'Accelerating'}, ...
            'ValueChangedFcn', @on_scenario_change);

        uilabel(g, 'Text', 'Noise');
        S.noiseDrop = uidropdown(g, ...
            'Items', {'None', '20 dB', '10 dB', '0 dB', '-5 dB'}, ...
            'ValueChangedFcn', @on_scenario_change);

        uilabel(g, 'Text', 'Target spacing');
        S.spacingDrop = uidropdown(g, 'Items', {'normal', 'close'}, ...
            'ValueChangedFcn', @on_scenario_change);

        uilabel(g, 'Text', 'Prediction');
        S.predictionLabel = uilabel(g, 'Text', '--', 'WordWrap', 'on');
    end

    function build_algorithm_panel(parent)

        p = uipanel(parent, 'Title', 'Image formation');
        g = uigridlayout(p, [3 1]);
        g.RowHeight = repmat({24}, 1, 3);

        S.ompCheck  = uicheckbox(g, 'Text', 'OMP',            'Value', true);
        S.nompCheck = uicheckbox(g, 'Text', 'NOMP',           'Value', true);
        S.bpCheck   = uicheckbox(g, 'Text', 'Backprojection', 'Value', true);
    end

    function build_advanced_panel(parent)

        p = uipanel(parent, 'Title', 'Radar, target and solver');
        g = uigridlayout(p, [16 2]);
        g.ColumnWidth = {170, '1x'};
        g.RowHeight   = repmat({26}, 1, 16);

        S.numFields = struct();

        add_num(g, 'Ks',                  'Number of scatterers',       10);
        add_num(g, 'oversampling_factor', 'Oversampling factor',         4);
        add_num(g, 'N_critical',          'Critical range pixels',      41);
        add_num(g, 'Nd',                  'Phase-history dimension',    16);
        add_num(g, 'prf',                 'PRF [Hz]',                  200);
        add_num(g, 'fc_GHz',              'Center frequency [GHz]',      1);
        add_num(g, 'fs_MHz',              'Sampling frequency [MHz]',  300);
        add_num(g, 'u0',                  'Range to rotation center [m]', 1000);
        add_num(g, 'w0',                  'Yaw rate w0 [rad/s]',        pi);
        add_num(g, 'w1',                  'Yaw accel w1 [rad/s^2]',    1e3);
        add_num(g, 'w2',                  'Yaw jerk w2 [rad/s^3]',     1e3);
        add_num(g, 'Rs',                  'NOMP Newton steps Rs',        4);
        add_num(g, 'Rc',                  'NOMP cyclic refinements Rc',  2);
        add_num(g, 'seed',                'RNG seed',                    0);

        uilabel(g, 'Text', 'Maneuvering trajectory');
        S.maneuverCheck = uicheckbox(g, 'Text', 'when accelerating', ...
            'Value', true, 'ValueChangedFcn', @(~,~) render_trajectory());

        uilabel(g, 'Text', 'Report rank(A)');
        S.rankCheck = uicheckbox(g, 'Text', 'slow for dense grids', 'Value', false);
    end

    function add_num(g, key, label, default)
        uilabel(g, 'Text', label);
        S.numFields.(key) = uieditfield(g, 'numeric', 'Value', default, ...
            'ValueChangedFcn', @(~,~) render_trajectory());
    end

    function build_run_controls(parent)

        p = uipanel(parent, 'Title', 'Run');
        g = uigridlayout(p, [4 2]);
        g.ColumnWidth = {'1x', '1x'};
        g.RowHeight   = {32, 26, 26, 'fit'};

        S.runButton = uibutton(g, 'Text', 'Run scenario', ...
            'FontWeight', 'bold', 'ButtonPushedFcn', @on_run);
        S.runButton.Layout.Column = [1 2];

        uilabel(g, 'Text', 'Display');
        S.logScaleCheck = uicheckbox(g, 'Text', 'log scale [dB]', 'Value', true);

        uilabel(g, 'Text', 'Dynamic range [dB]');
        S.dynRangeField = uieditfield(g, 'numeric', 'Value', 60, ...
            'ValueChangedFcn', @(~,~) redraw_images());

        S.popButton = uibutton(g, 'Text', 'Pop out figure', ...
            'ButtonPushedFcn', @on_popout, 'Enable', 'off');
        S.exportButton = uibutton(g, 'Text', 'Export history', ...
            'ButtonPushedFcn', @on_export, 'Enable', 'off');

        S.logScaleCheck.ValueChangedFcn = @(~,~) redraw_images();
    end

    function build_trajectory_tab()

        g = uigridlayout(S.trajTab, [2 1]);
        g.RowHeight   = {'1x', '1x'};
        g.ColumnWidth = {'1x'};

        % the derived quantities ride along as axis subtitles rather than as
        % a separate label: a component built inside a tab that is not the
        % selected one does not get positioned by the grid
        S.thetaAxes = uiaxes(g);
        S.thetaAxes.Layout.Row = 1;

        S.rateAxes = uiaxes(g);
        S.rateAxes.Layout.Row = 2;
    end

    function build_error_tab()

        g = uigridlayout(S.errorTab, [4 1]);
        g.RowHeight = {'fit', 'fit', 'fit', '1x'};

        uilabel(g, 'Text', 'Per-algorithm summary', 'FontWeight', 'bold');
        S.summaryTable = uitable(g, 'ColumnName', ...
            {'Algorithm', 'RMS error [m]', 'Matched', 'Missed', ...
             'False alarms', 'Median error [m]', 'Max error [m]'});
        S.summaryTable.Layout.Row = 2;

        uilabel(g, 'Text', 'Matched scatterer pairs', 'FontWeight', 'bold');
        S.pairTable = uitable(g, 'ColumnName', ...
            {'Algorithm', 'True crossrange [m]', 'True range [m]', ...
             'Est crossrange [m]', 'Est range [m]', 'Error [m]'});
    end

    function build_history_tab()

        g = uigridlayout(S.historyTab, [2 1]);
        g.RowHeight = {'fit', '1x'};

        uilabel(g, 'Text', ...
            ['Every run this session. Change one dropdown, run again, and ' ...
             'compare the rows.'], 'WordWrap', 'on');

        S.historyTable = uitable(g, 'ColumnName', history_columns());
    end

    function placeholder(g, text)
        delete(g.Children);
        g.RowHeight    = {'1x'};
        g.ColumnWidth  = {'1x'};
        lbl = uilabel(g, 'Text', text, 'HorizontalAlignment', 'center', ...
            'FontColor', [0.45 0.45 0.45]);
        lbl.Layout.Row = 1;
        lbl.Layout.Column = 1;
    end

%% ------------------------------------------------------- preset wiring

    function apply_preset_items()

        if isempty(S.scenarios)
            S.presetDrop.Items = {'(custom)'};
            S.presetDrop.ItemsData = 0;
            return
        end

        t = S.scenarios;
        items = cell(1, height(t) + 1);
        items{1} = '(custom)';
        for i = 1:height(t)
            items{i+1} = sprintf('Case %d  --  %d/%d %s %s %s %s', ...
                t.Case_(i), ...
                t.NumberOfAmbiguitiesHavingScatterers(i), ...
                t.NumberOfAmbiguitiesInImageFormer(i), ...
                char(string(t.ScattererLocations(i))), ...
                char(string(t.ImageFormerGridDensity(i))), ...
                char(string(t.AngleRate(i))), ...
                char(string(t.TargetSpacing(i))));
        end

        S.presetDrop.Items     = items;
        S.presetDrop.ItemsData = 0:height(t);
        S.presetDrop.Value     = 0;
    end

    function on_preset(~, ~)

        row = S.presetDrop.Value;
        if row == 0 || isempty(S.scenarios)
            S.predictionLabel.Text = '--';
            return
        end

        t = S.scenarios(row, :);

        S.ambScatDrop.Value = t.NumberOfAmbiguitiesHavingScatterers;
        S.ambIfDrop.Value   = t.NumberOfAmbiguitiesInImageFormer;

        % set_choice widens a dropdown rather than erroring, so a value added
        % to the spreadsheet later (another Noise level, say) still loads
        set_choice(S.scatLocDrop,     canonical_grid_label(t.ScattererLocations));
        set_choice(S.gridDensityDrop, t.ImageFormerGridDensity);
        set_choice(S.angleRateDrop,   t.AngleRate);
        set_choice(S.noiseDrop,       t.Noise);
        set_choice(S.spacingDrop,     t.TargetSpacing);

        S.predictionLabel.Text  = char(string(t.Prediction));
        sync_enable();
        render_trajectory();
    end

    function on_scenario_change(~, ~)
        % any hand edit means we are no longer on a numbered case
        S.presetDrop.Value = 0;
        S.predictionLabel.Text = '--';
        sync_enable();
        render_trajectory();
    end

    function sync_enable()
        % the oversampling factor is only meaningful on an oversampled grid
        is_over = strcmpi(S.gridDensityDrop.Value, 'Oversampled');
        S.numFields.oversampling_factor.Enable = matlab.lang.OnOffSwitchState(is_over);

        % closely spaced targets are only synthesized for off-grid scatterers
        is_off = strcmpi(S.scatLocDrop.Value, 'Off-Grid');
        S.spacingDrop.Enable = matlab.lang.OnOffSwitchState(is_off);

        is_accel = strcmpi(S.angleRateDrop.Value, 'Accelerating');
        S.maneuverCheck.Enable = matlab.lang.OnOffSwitchState(is_accel);
        S.numFields.w1.Enable  = matlab.lang.OnOffSwitchState(is_accel);
        S.numFields.w2.Enable  = matlab.lang.OnOffSwitchState(is_accel);
    end

%% ---------------------------------------------------------------- run

    function on_run(~, ~)

        cfg  = collect_cfg();
        opts = collect_opts();

        if ~(opts.execute_omp || opts.execute_nomp || opts.execute_bp)
            uialert(S.fig, 'Select at least one image formation algorithm.', ...
                'Nothing to run');
            return
        end

        S.logArea.Value = {''};
        S.runButton.Enable = 'off';
        dlg = uiprogressdlg(S.fig, 'Title', 'Running scenario', ...
            'Message', 'starting', 'Indeterminate', 'off', 'Value', 0);
        restore = onCleanup(@() cleanup_run(dlg));

        opts.log_fcn      = @append_log;
        opts.progress_fcn = @(frac, msg) update_progress(dlg, frac, msg);

        try
            out = isar_run_scenario(cfg, opts);
        catch ME
            append_log(['ERROR: ' ME.message]);
            uialert(S.fig, ME.message, 'Scenario failed');
            return
        end

        S.last = out;
        render_trajectory();
        redraw_images();
        render_errors(out);
        append_history(out);

        S.popButton.Enable    = 'on';
        S.exportButton.Enable = 'on';
    end

    function cleanup_run(dlg)
        if isvalid(dlg)
            close(dlg);
        end
        S.runButton.Enable = 'on';
    end

    function update_progress(dlg, frac, msg)
        if ~isvalid(dlg)
            return
        end
        dlg.Value   = max(0, min(1, frac));
        dlg.Message = msg;
        drawnow limitrate
    end

    function append_log(text)
        lines = S.logArea.Value;
        if numel(lines) == 1 && isempty(lines{1})
            lines = {};
        end
        S.logArea.Value = [lines(:); {text}];
        drawnow limitrate
    end

    function cfg = collect_cfg()
        cfg = struct( ...
            'NumberOfAmbiguitiesHavingScatterers', S.ambScatDrop.Value, ...
            'NumberOfAmbiguitiesInImageFormer',    S.ambIfDrop.Value, ...
            'ScattererLocations',                  S.scatLocDrop.Value, ...
            'ImageFormerGridDensity',              S.gridDensityDrop.Value, ...
            'AngleRate',                           S.angleRateDrop.Value, ...
            'Noise',                               S.noiseDrop.Value, ...
            'TargetSpacing',                       S.spacingDrop.Value);
    end

    function opts = collect_opts()
        n = S.numFields;
        opts = struct( ...
            'execute_omp',         S.ompCheck.Value, ...
            'execute_nomp',        S.nompCheck.Value, ...
            'execute_bp',          S.bpCheck.Value, ...
            'Ks',                  n.Ks.Value, ...
            'oversampling_factor', n.oversampling_factor.Value, ...
            'N_critical',          n.N_critical.Value, ...
            'Nd',                  n.Nd.Value, ...
            'prf',                 n.prf.Value, ...
            'fc',                  n.fc_GHz.Value * 1e9, ...
            'fs',                  n.fs_MHz.Value * 1e6, ...
            'u0',                  n.u0.Value, ...
            'w0',                  n.w0.Value, ...
            'w1',                  n.w1.Value, ...
            'w2',                  n.w2.Value, ...
            'Rs',                  n.Rs.Value, ...
            'Rc',                  n.Rc.Value, ...
            'seed',                n.seed.Value, ...
            'complex_maneuver',    S.maneuverCheck.Value, ...
            'compute_rank',        S.rankCheck.Value);
    end

%% ------------------------------------------------------------ plotting

    function redraw_images()

        if isempty(S.last)
            return
        end

        out    = S.last;
        panels = result_panels(out.x_hat);
        n      = size(panels, 1);

        delete(S.imageGrid.Children);
        [rows, cols] = panel_geometry(n);
        S.imageGrid.RowHeight   = repmat({'1x'}, 1, rows);
        S.imageGrid.ColumnWidth = repmat({'1x'}, 1, cols);

        S.axesList = gobjects(1, n);
        for i = 1:n
            ax = uiaxes(S.imageGrid);
            ax.Layout.Row    = ceil(i / cols);
            ax.Layout.Column = mod(i - 1, cols) + 1;
            draw_panel(ax, panels{i,1}, panels{i,2}, out, ...
                S.logScaleCheck.Value, S.dynRangeField.Value);
            S.axesList(i) = ax;
        end
    end

    function render_trajectory()
    % The trajectory depends only on the scenario and the motion settings,
    % not on the reconstruction, so it is redrawn live as those change --
    % you can see the motion before paying for a run.

        if ~isfield(S, 'thetaAxes') || ~isvalid(S.thetaAxes)
            return
        end

        try
            traj = isar_target_trajectory(collect_cfg(), collect_opts());
        catch ME
            trajectory_message(['trajectory unavailable: ' ME.message]);
            return
        end

        if isempty(traj.t_m)
            trajectory_message( ...
                'no pulses -- check the phase-history dimension and PRF');
            return
        end

        t_ms = traj.t_m * 1e3;

        % markers sit on the pulses themselves, so the slow-time sampling is
        % as visible as the underlying motion
        plot(S.thetaAxes, t_ms, traj.theta, '-o', 'LineWidth', 1.5, ...
            'MarkerSize', 4, 'MarkerFaceColor', [0 0.447 0.741]);
        grid(S.thetaAxes, 'on');
        xlabel(S.thetaAxes, 'slow time [ms]');
        ylabel(S.thetaAxes, '\theta [rad]');
        title(S.thetaAxes, 'Target yaw angle vs slow time');

        plot(S.rateAxes, t_ms, traj.theta_rate, '-o', 'LineWidth', 1.5, ...
            'MarkerSize', 4, 'MarkerFaceColor', [0.85 0.325 0.098], ...
            'Color', [0.85 0.325 0.098]);
        grid(S.rateAxes, 'on');
        xlabel(S.rateAxes, 'slow time [ms]');
        ylabel(S.rateAxes, 'd\theta/dt [rad/s]');
        title(S.rateAxes, 'Target yaw rate vs slow time');

        if numel(t_ms) > 1
            xlim(S.thetaAxes, [t_ms(1) t_ms(end)]);
            xlim(S.rateAxes,  [t_ms(1) t_ms(end)]);
        end

        span     = max(traj.theta) - min(traj.theta);
        sin_span = max(sin(traj.theta)) - min(sin(traj.theta));

        if traj.maneuvering
            model = 'piecewise maneuvering trajectory (sign-alternating jerk)';
        else
            model = 'polynomial yaw: theta = w0 t + w1 t^2/2 + w2 t^3/3';
        end

        % the angular span sets the cross-range resolution, so it is worth
        % reading off the same panel as the motion that produced it
        subtitle(S.thetaAxes, sprintf( ...
            '%d pulses over %.2f ms  |  span %.4f rad (%.2f deg)  |  sin span %.4f', ...
            numel(traj.t_m), traj.T*1e3, span, rad2deg(span), sin_span), ...
            'Interpreter', 'none');

        subtitle(S.rateAxes, sprintf('%.3f to %.3f rad/s  |  %s', ...
            min(traj.theta_rate), max(traj.theta_rate), model), ...
            'Interpreter', 'none');
    end

    function trajectory_message(msg)
        cla(S.thetaAxes);
        cla(S.rateAxes);
        title(S.thetaAxes, msg, 'Interpreter', 'none');
        subtitle(S.thetaAxes, '');
        title(S.rateAxes, '');
        subtitle(S.rateAxes, '');
    end

    function on_popout(~, ~)

        if isempty(S.last)
            return
        end

        out    = S.last;
        panels = result_panels(out.x_hat);
        n      = size(panels, 1);

        [rows, cols] = panel_geometry(n);
        f = figure('Name', scenario_title(out.cfg), ...
            'Position', [100 100 520*cols 520*rows]);
        for i = 1:n
            ax = subplot(rows, cols, i, 'Parent', f);
            draw_panel(ax, panels{i,1}, panels{i,2}, out, ...
                S.logScaleCheck.Value, S.dynRangeField.Value);
        end
        sgtitle(f, scenario_title(out.cfg), 'Interpreter', 'none');
    end

%% ----------------------------------------------------------- results

    function render_errors(out)

        algs = present_algorithms(out.x_hat);

        summary = cell(numel(algs), 7);
        pairs   = {};

        for i = 1:numel(algs)
            alg = algs{i};
            r   = out.x_hat.(alg);

            summary(i, :) = { upper(alg), r.error, size(r.pairs, 1), ...
                numel(r.missed), numel(r.false_alarms), ...
                median(r.d), max(r.d) };

            for j = 1:size(r.pairs, 1)
                it = r.pairs(j, 1);   % index into the true locations
                ie = r.pairs(j, 2);   % index into the estimated locations
                pairs(end+1, :) = { upper(alg), ...
                    out.latent_locations(it, 1), out.latent_locations(it, 2), ...
                    r.positions(ie, 1),          r.positions(ie, 2), ...
                    r.d(j) };  %#ok<AGROW>
            end
        end

        S.summaryTable.Data = summary;
        S.pairTable.Data    = pairs;
    end

    function append_history(out)

        algs = present_algorithms(out.x_hat);
        err  = @(a) pick_error(out.x_hat, a, algs);

        % errors lead so a Critical vs Oversampled comparison is readable
        % without scrolling the table sideways
        row = { size(S.history, 1) + 1, ...
            err('omp'), err('nomp'), err('bp'), ...
            out.cfg.NumberOfAmbiguitiesHavingScatterers, ...
            out.cfg.NumberOfAmbiguitiesInImageFormer, ...
            out.cfg.ScattererLocations, ...
            out.cfg.ImageFormerGridDensity, ...
            out.cfg.AngleRate, ...
            out.cfg.Noise, ...
            out.cfg.TargetSpacing, ...
            out.sim_config.Nx, ...
            out.sim_config.Ny, ...
            out.sim_config.CrossRangePixelRes_m, ...
            out.sim_config.RangePixelRes_m, ...
            char(out.sim_config.DopplerAliasing) };

        S.history(end+1, :) = row;
        S.historyTable.Data = S.history;
    end

    function on_export(~, ~)

        if isempty(S.history)
            return
        end

        [file, path] = uiputfile('*.csv', 'Export run history', ...
            fullfile(S.root, 'results', 'gui_history.csv'));
        if isequal(file, 0)
            return
        end

        t = cell2table(S.history, 'VariableNames', ...
            matlab.lang.makeValidName(history_columns()));
        writetable(t, fullfile(path, file));
        append_log(sprintf('exported %d runs to %s', ...
            size(S.history, 1), fullfile(path, file)));
    end
end

%% ================================================== local functions ==

function scenarios = load_scenarios(root)
% LOAD_SCENARIOS  Read the scenario matrix, tolerating its absence.

    scenarios = [];
    file = fullfile(root, 'documentation', 'scenarios_new.xlsx');
    if ~isfile(file)
        return
    end

    try
        scenarios = readtable(file);
    catch
        scenarios = [];
    end
end

function set_choice(dropdown, value)
% SET_CHOICE  Select VALUE on DROPDOWN, adding it to the list if it is new.

    value = char(string(value));

    match = strcmpi(dropdown.Items, value);
    if any(match)
        dropdown.Value = dropdown.Items{find(match, 1)};
    else
        dropdown.Items = [dropdown.Items, {value}];
        dropdown.Value = value;
    end
end

function label = canonical_grid_label(value)
% CANONICAL_GRID_LABEL  The spreadsheet mixes 'Off-Grid' and 'Off-grid'.

    if strcmpi(char(string(value)), 'off-grid')
        label = 'Off-Grid';
    else
        label = 'On-Grid';
    end
end

function algs = present_algorithms(x_hat)
% PRESENT_ALGORITHMS  Which algorithms produced a scored result.

    algs = {};
    for a = {'omp', 'nomp', 'bp'}
        if is_scored(x_hat, a{1})
            algs{end+1} = a{1}; %#ok<AGROW>
        end
    end
end

function tf = is_scored(x_hat, alg)
    tf = isfield(x_hat, alg) && isfield(x_hat.(alg), 'error');
end

function panels = result_panels(x_hat)
% RESULT_PANELS  The plots to show, as {algorithm, mode} rows.
%
%   An algorithm can contribute more than one view. OMP gets both: the
%   gridded image, and the extracted peak positions drawn the same way as
%   NOMP's, which is the easier of the two to compare against truth when
%   the scatterers are closely spaced.

    panels = cell(0, 2);

    if is_scored(x_hat, 'omp')
        panels(end+1, :) = {'omp', 'image'};
        panels(end+1, :) = {'omp', 'positions'};
    end

    if is_scored(x_hat, 'nomp')
        panels(end+1, :) = {'nomp', 'positions'};
    end

    if is_scored(x_hat, 'bp')
        panels(end+1, :) = {'bp', 'image'};
    end
end

function [rows, cols] = panel_geometry(n)
    rows = 1 + (n > 2);
    cols = ceil(n / max(rows, 1));
end

function e = pick_error(x_hat, alg, algs)
    if ismember(alg, algs)
        e = x_hat.(alg).error;
    else
        e = NaN;
    end
end

function cols = history_columns()
    cols = {'Run', 'RMS_OMP_m', 'RMS_NOMP_m', 'RMS_BP_m', ...
        'AmbWithScatterers', 'AmbInImageFormer', 'ScattererLocations', ...
        'GridDensity', 'AngleRate', 'Noise', 'TargetSpacing', ...
        'Nx', 'Ny', 'CrossRangePixelRes_m', 'RangePixelRes_m', ...
        'DopplerAliasing'};
end

function t = scenario_title(cfg)
    t = sprintf('%d/%d amb, %s, %s, %s, %s spacing, noise %s', ...
        cfg.NumberOfAmbiguitiesHavingScatterers, ...
        cfg.NumberOfAmbiguitiesInImageFormer, ...
        cfg.ScattererLocations, cfg.ImageFormerGridDensity, ...
        cfg.AngleRate, cfg.TargetSpacing, cfg.Noise);
end

function draw_panel(ax, alg, mode, out, log_scale, dyn_range)
% DRAW_PANEL  Render one view of one algorithm's result into AX.
%
%   mode 'image'     -- the gridded latent image, truth overlaid
%   mode 'positions' -- the recovered scatterer positions alone, truth
%                       overlaid, identically for OMP and NOMP
%
%   Works with both uiaxes (in the GUI) and ordinary axes (in a popped out
%   figure), so the two always show the same thing.

    x_array = out.x_array;
    y_array = out.y_array;
    u0      = out.u0;
    r       = out.x_hat.(alg);

    true_x = out.latent_locations(:, 1);
    true_y = out.latent_locations(:, 2) + u0;

    cla(ax);

    switch mode
        case 'image'
            img = abs(r.image);

            if log_scale
                img = 20*log10(img);
                top = max(img(isfinite(img)));
                if isempty(top)
                    top = 0;
                end
                imagesc(ax, x_array, y_array + u0, img);
                clim(ax, [top - dyn_range, top]);
                cbar_label = 'Log-Scaled [dB]';
            else
                imagesc(ax, x_array, y_array + u0, img);
                cbar_label = 'Amplitude [Linear]';
            end

            colormap(ax, gray);
            cb = colorbar(ax);
            cb.Label.String = cbar_label;
            heading = upper(alg);

        case 'positions'
            % OMP's positions are the extracted image peaks and NOMP's are
            % estimated off-grid directly, but both are [x, y] in metres, so
            % they are drawn the same way and read the same way
            plot(ax, r.positions(:, 1), r.positions(:, 2) + u0, '*', ...
                'MarkerEdgeColor', [0 0 1], 'MarkerSize', 10, ...
                'LineWidth', 1.5);
            set(ax, 'YDir', 'reverse');
            heading = [upper(alg) ' positions'];

        otherwise
            error('isar_gui:panelMode', 'unknown panel mode ''%s''', mode);
    end

    hold(ax, 'on');
    plot(ax, true_x, true_y, 'o', 'MarkerEdgeColor', [1 0 0], ...
        'MarkerSize', 10, 'LineWidth', 1.5);
    hold(ax, 'off');

    axis(ax, 'square');
    xlim(ax, [min(x_array), max(x_array)]);
    ylim(ax, [min(y_array) + u0, max(y_array) + u0]);
    xlabel(ax, 'crossrange [m]');
    ylabel(ax, 'range [m]');
    title(ax, sprintf('%s  (RMS error: %.3f m)', heading, r.error));
end
