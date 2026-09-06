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
%   Run it as  gui/isar_gui  from the project root, or add gui/ to the path.
%   It resolves documentation/, results/ and functions/ relative to its own
%   location, so the working directory does not matter.
%
%   Every run is appended to the History tab, so stepping a single dropdown
%   -- Critical to Oversampled, say -- builds a side-by-side comparison you
%   can export to CSV.

    % this file lives in gui/, so the project root is one level up. both the
    % root and functions/ go on the path: the root holds the scenario helper
    % scripts, functions/ holds everything the solvers need
    root = fileparts(fileparts(mfilename('fullpath')));
    addpath(root);
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

        % the gridded images and the extracted positions answer different
        % questions and are read differently, so they get a tab each rather
        % than competing for room in one grid
        S.positionTab = uitab(S.tabs, 'Title', 'Scatterer positions');
        S.positionGrid = uigridlayout(S.positionTab, [1 1]);
        placeholder(S.positionGrid, 'Press Run to recover positions.');

        S.trajTab = uitab(S.tabs, 'Title', 'Trajectory');
        build_trajectory_tab();

        S.errorTab = uitab(S.tabs, 'Title', 'Position errors');
        build_error_tab();

        S.pathTab = uitab(S.tabs, 'Title', 'Refinement path');
        build_path_tab();

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
            'ItemsData', [1 3], 'Value', 3, ...
            'ValueChangedFcn', @on_scenario_change);

        uilabel(g, 'Text', '# ambiguities in image former');
        S.ambIfDrop = uidropdown(g, 'Items', {'1', '3'}, ...
            'ItemsData', [1 3], 'ValueChangedFcn', @on_scenario_change);

        uilabel(g, 'Text', 'Scatterer locations');
        S.scatLocDrop = uidropdown(g, 'Items', {'On-Grid', 'Off-Grid'}, ...
            'ValueChangedFcn', @on_scenario_change);

        uilabel(g, 'Text', 'Image former grid density');
        S.gridDensityDrop = uidropdown(g, 'Items', {'Critical', 'Oversampled'}, ...
            'Value', 'Oversampled', 'ValueChangedFcn', @on_scenario_change);

        uilabel(g, 'Text', 'Angle rate');
        S.angleRateDrop = uidropdown(g, 'Items', {'Constant', 'Accelerating'}, ...
            'Value', 'Accelerating', 'ValueChangedFcn', @on_scenario_change);

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
        g = uigridlayout(p, [5 1]);
        g.RowHeight = repmat({24}, 1, 5);

        S.ompCheck    = uicheckbox(g, 'Text', 'OMP',            'Value', true);
        S.modOmpCheck = uicheckbox(g, 'Text', 'Modified OMP',   'Value', true);
        S.nompCheck   = uicheckbox(g, 'Text', 'NOMP',           'Value', true);
        S.prompCheck  = uicheckbox(g, 'Text', 'PROMP',          'Value', true);
        S.bpCheck     = uicheckbox(g, 'Text', 'Backprojection', 'Value', true);
    end

    function build_advanced_panel(parent)

        p = uipanel(parent, 'Title', 'Radar, target and solver');
        g = uigridlayout(p, [21 2]);
        g.ColumnWidth = {170, '1x'};
        g.RowHeight   = repmat({26}, 1, 21);

        S.numFields = struct();

        add_num(g, 'Ks',                  'Number of scatterers',        2);
        add_num(g, 'oversampling_factor', 'Oversampling factor',         4);
        add_num(g, 'N_critical',          'Critical range pixels',      15);
        add_num(g, 'Nd',                  'Phase-history dimension',    16);
        add_num(g, 'prf',                 'PRF [Hz]',                 6000);
        add_num(g, 'fc_GHz',              'Center frequency [GHz]',     30);
        add_num(g, 'fs_MHz',              'Sampling frequency [MHz]',  300);
        add_num(g, 'u0',                  'Range to rotation center [m]', 1000);
        add_num(g, 'w0',                  'Yaw rate w0 [rad/s]',        pi);
        add_num(g, 'w1',                  'Yaw accel w1 [rad/s^2]',    170);
        add_num(g, 'w2',                  'Yaw jerk w2 [rad/s^3]',       0);
        add_num(g, 'Rs',                  'Newton/GN steps Rs',          4);
        add_num(g, 'Rc',                  'NOMP cyclic refinements Rc',  2);
        add_num(g, 'amb_refine_pixels',   'Mod-OMP search [pixels]',     5);
        add_num(g, 'seed',                'RNG seed',                    0);

        uilabel(g, 'Text', 'Maneuvering trajectory');
        S.maneuverCheck = uicheckbox(g, 'Text', 'when accelerating', ...
            'Value', false, 'ValueChangedFcn', @(~,~) render_trajectory());

        uilabel(g, 'Text', 'Report rank(A)');
        S.rankCheck = uicheckbox(g, 'Text', 'slow for dense grids', 'Value', false);

        uilabel(g, 'Text', 'Range cell only');
        S.rangeCellCheck = uicheckbox(g, 'Text', '1-D crossrange', ...
            'Value', false, 'ValueChangedFcn', @(~,~) sync_enable());

        add_num(g, 'range_cell',          'Range cell [m]',              0);

        uilabel(g, 'Text', 'Refinement path');
        S.histCheck = uicheckbox(g, 'Text', 'record refinement steps', 'Value', true);

        add_num(g, 'save_atom_idx',       'Traced atom index',           1);
    end

    function add_num(g, key, label, default)
        uilabel(g, 'Text', label);
        S.numFields.(key) = uieditfield(g, 'numeric', 'Value', default, ...
            'ValueChangedFcn', @(~,~) render_trajectory());
    end

    function build_run_controls(parent)

        p = uipanel(parent, 'Title', 'Run');
        g = uigridlayout(p, [10 2]);
        g.ColumnWidth = {'1x', '1x'};
        g.RowHeight   = {32, 26, 26, 26, 26, 26, 26, 26, 26, 'fit'};

        S.runButton = uibutton(g, 'Text', 'Run scenario', ...
            'FontWeight', 'bold', 'ButtonPushedFcn', @on_run);
        S.runButton.Layout.Column = [1 2];

        uilabel(g, 'Text', 'Display');
        S.logScaleCheck = uicheckbox(g, 'Text', 'log scale [dB]', 'Value', true);

        uilabel(g, 'Text', 'Dynamic range [dB]');
        S.dynRangeField = uieditfield(g, 'numeric', 'Value', 60, ...
            'ValueChangedFcn', @(~,~) redraw_images());

        uilabel(g, 'Text', 'Ambiguities');
        S.ambViewCheck = uicheckbox(g, 'Text', 'show -1, 0, +1', ...
            'Value', false, 'ValueChangedFcn', @(~,~) redraw_images());

        uilabel(g, 'Text', 'Aliased targets');
        S.aliasCheck = uicheckbox(g, 'Text', 'fold into ambiguity 0', ...
            'Value', true, 'ValueChangedFcn', @(~,~) redraw_images());

        uilabel(g, 'Text', 'Save output');
        S.saveCheck = uicheckbox(g, 'Text', 'after each run', 'Value', false);

        uilabel(g, 'Text', 'Folder');
        S.saveFolderField = uieditfield(g, 'text', ...
            'Value', fullfile(S.root, 'results'));

        S.browseButton = uibutton(g, 'Text', 'Browse...', ...
            'ButtonPushedFcn', @on_browse_save);
        S.saveNowButton = uibutton(g, 'Text', 'Save last run', ...
            'ButtonPushedFcn', @on_save_now, 'Enable', 'off');

        uilabel(g, 'Text', 'Include');
        S.saveDictsCheck = uicheckbox(g, 'Text', 'dictionaries (2nd file)', ...
            'Value', false);

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

    function build_path_tab()
    % The trace of NOMP's Newton steps for one atom, zoomed on the true
    % scatterer it converged to, alongside how fast it got there.

        g = uigridlayout(S.pathTab, [2 2]);
        g.ColumnWidth = {'2x', '1x'};
        g.RowHeight   = {'1x', '1x'};

        % one row per algorithm, trace on the left and convergence on the
        % right, so NOMP and PROMP can be read against each other directly
        S.pathAlgs = {'nomp', 'promp'};
        S.pathAxes = gobjects(1, 2);
        S.convAxes = gobjects(1, 2);
        for ip = 1:2
            S.pathAxes(ip) = uiaxes(g);
            S.pathAxes(ip).Layout.Row = ip;  S.pathAxes(ip).Layout.Column = 1;
            S.convAxes(ip) = uiaxes(g);
            S.convAxes(ip).Layout.Row = ip;  S.convAxes(ip).Layout.Column = 2;
        end
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

        % a single range cell has no range axis, so the range-grid settings
        % and the range-density choice have nothing to act on
        one_cell = S.rangeCellCheck.Value;
        S.numFields.N_critical.Enable = matlab.lang.OnOffSwitchState(~one_cell);
        S.numFields.range_cell.Enable = matlab.lang.OnOffSwitchState(one_cell);

        is_accel = strcmpi(S.angleRateDrop.Value, 'Accelerating');
        S.maneuverCheck.Enable = matlab.lang.OnOffSwitchState(is_accel);
        S.numFields.w1.Enable  = matlab.lang.OnOffSwitchState(is_accel);
        S.numFields.w2.Enable  = matlab.lang.OnOffSwitchState(is_accel);
    end

%% ---------------------------------------------------------------- run

    function on_run(~, ~)

        cfg  = collect_cfg();
        opts = collect_opts();

        if ~(opts.execute_omp || opts.execute_nomp || opts.execute_promp ...
                || opts.execute_mod_omp || opts.execute_bp)
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
        redraw_refinement_path();
        render_errors(out);
        append_history(out);

        S.popButton.Enable     = 'on';
        S.exportButton.Enable  = 'on';
        S.saveNowButton.Enable = 'on';

        if S.saveCheck.Value
            save_run_output(out);
        end
    end

    function on_browse_save(~, ~)

        start = strtrim(S.saveFolderField.Value);
        if isempty(start) || ~isfolder(start)
            start = S.root;
        end

        folder = uigetdir(start, 'Choose a folder for saved runs');
        if isequal(folder, 0)
            return
        end
        S.saveFolderField.Value = folder;
    end

    function on_save_now(~, ~)

        if isempty(S.last)
            uialert(S.fig, 'Run a scenario first.', 'Nothing to save');
            return
        end
        save_run_output(S.last);
    end

    function save_run_output(out)
    % SAVE_RUN_OUTPUT  Write one run to <folder>/isar_run_<timestamp>.{mat,png}
    %
    %   The .mat holds one struct, run_data, carrying x_hat -- the OMP and BP
    %   images and the NOMP and PROMP positions and reflectivities -- together
    %   with the grids, the truth, the measurement, and the configuration. The
    %   .png is the same panel figure that 'Pop out figure' produces.
    %
    %   With 'dictionaries' ticked, a second file <stem>_dictionaries.mat is
    %   written holding A, as and y as top-level variables, so that load() puts
    %   them straight into the workspace. They live apart from the run file
    %   because they are one to three orders of magnitude larger than
    %   everything else put together.

        folder = strtrim(S.saveFolderField.Value);
        if isempty(folder)
            folder = fullfile(S.root, 'results');
        end

        if ~isfolder(folder)
            [ok, msg] = mkdir(folder);
            if ~ok
                uialert(S.fig, msg, 'Could not create folder');
                return
            end
        end

        stamp = char(string(datetime('now', 'Format', 'yyyyMMdd_HHmmss')));
        stem  = fullfile(folder, ['isar_run_' stamp]);

        run_data                  = struct();
        run_data.cfg              = out.cfg;
        run_data.sim_config       = out.sim_config;
        run_data.x_hat            = out.x_hat;
        run_data.x_array          = out.x_array;
        run_data.y_array          = out.y_array;
        run_data.u0               = out.u0;
        run_data.theta            = out.theta;
        run_data.t_m              = out.t_m;
        run_data.target_locations = out.target_locations;
        run_data.latent_locations = out.latent_locations;
        run_data.amb_of_k         = out.amb_of_k;
        run_data.is_latent        = out.is_latent;

        if isfield(out, 'y')
            run_data.measurement = out.y;
        end

        % opts carries handles to nested functions, which cannot be saved
        run_data.opts = rmfield(out.opts, intersect(fieldnames(out.opts), ...
            {'log_fcn', 'progress_fcn'}));

        try
            save([stem '.mat'], 'run_data');
        catch ME
            append_log(['ERROR saving .mat: ' ME.message]);
            uialert(S.fig, ME.message, 'Save failed');
            return
        end

        % the dictionaries dwarf everything else, so they go to their own file
        % and as top-level variables -- load(f) puts A, as and y straight into
        % the workspace, and the run file above stays small enough to reload
        % many times when comparing runs
        dict_file = '';
        if S.saveDictsCheck.Value
            A  = out.A;    % reconstruction dictionary  [ML x Nx*Ny]
            as = out.as;   % measurement sensing matrix [ML x Ks]
            y  = out.y;    % measurement                [ML x 1]

            dict_file = [stem '_dictionaries.mat'];
            sz = whos('A', 'as', 'y');
            try
                % MAT-file v7 tops out at 2 GB
                if sum([sz.bytes]) > 1.9e9
                    save(dict_file, 'A', 'as', 'y', '-v7.3');
                else
                    save(dict_file, 'A', 'as', 'y');
                end
            catch ME
                append_log(['ERROR saving dictionaries: ' ME.message]);
                uialert(S.fig, ME.message, 'Save failed');
                return
            end
        end

        % the formed images, rendered the same way as the popped out figure
        panels = result_panels(out.x_hat);
        n      = size(panels, 1);
        [rows, cols] = panel_geometry(n);

        f = figure('Visible', 'off', 'Position', [100 100 520*cols 520*rows]);
        restore_fig = onCleanup(@() close(f));
        for i = 1:n
            ax = subplot(rows, cols, i, 'Parent', f);
            draw_panel(ax, panels{i,1}, panels{i,2}, out, ...
                S.logScaleCheck.Value, S.dynRangeField.Value, ...
                S.ambViewCheck.Value, S.aliasCheck.Value);
        end
        sgtitle(f, scenario_title(out.cfg), 'Interpreter', 'none');
        exportgraphics(f, [stem '.png'], 'Resolution', 150);

        info = dir([stem '.mat']);
        append_log(sprintf('saved %s (%.1f MB) and %s', ...
            [stem '.mat'], info.bytes/2^20, [stem '.png']));
        if ~isempty(dict_file)
            info = dir(dict_file);
            append_log(sprintf('saved %s (%.1f MB) holding A, as and y', ...
                dict_file, info.bytes/2^20));
        end
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
            'execute_promp',       S.prompCheck.Value, ...
            'execute_mod_omp',     S.modOmpCheck.Value, ...
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
            'amb_refine_pixels',   n.amb_refine_pixels.Value, ...
            'seed',                n.seed.Value, ...
            'complex_maneuver',    S.maneuverCheck.Value, ...
            'compute_rank',        S.rankCheck.Value, ...
            'range_cell_mode',     S.rangeCellCheck.Value, ...
            'range_cell',          n.range_cell.Value, ...
            'save_histories',      S.histCheck.Value, ...
            'save_atom_idx',       max(1, round(n.save_atom_idx.Value)));
    end

%% ------------------------------------------------------------ plotting

    function redraw_images()

        if isempty(S.last)
            return
        end

        S.axesList = [ ...
            render_panel_grid(S.imageGrid,    'image', ...
                'No gridded image was formed.'), ...
            render_panel_grid(S.positionGrid, 'positions', ...
                'No positions were recovered.')];
    end

    function axes_made = render_panel_grid(grid, mode, empty_text)
    % RENDER_PANEL_GRID  Draw every panel of one mode into one tab's grid.

        out    = S.last;
        panels = result_panels(out.x_hat, mode);
        n      = size(panels, 1);

        delete(grid.Children);

        if n == 0
            grid.RowHeight   = {'1x'};
            grid.ColumnWidth = {'1x'};
            placeholder(grid, empty_text);
            axes_made = gobjects(1, 0);
            return
        end

        [rows, cols] = panel_geometry(n);
        grid.RowHeight   = repmat({'1x'}, 1, rows);
        grid.ColumnWidth = repmat({'1x'}, 1, cols);

        axes_made = gobjects(1, n);
        for i = 1:n
            ax = uiaxes(grid);
            ax.Layout.Row    = ceil(i / cols);
            ax.Layout.Column = mod(i - 1, cols) + 1;
            draw_panel(ax, panels{i,1}, panels{i,2}, out, ...
                S.logScaleCheck.Value, S.dynRangeField.Value, ...
                S.ambViewCheck.Value, S.aliasCheck.Value);
            axes_made(i) = ax;
        end
    end

    function redraw_refinement_path()

        if isempty(S.last)
            return
        end

        drawn = false(1, numel(S.pathAlgs));
        ylims = zeros(0, 2);

        for ip = 1:numel(S.pathAlgs)
            alg = S.pathAlgs{ip};
            [drawn(ip), yl] = draw_refinement_path( ...
                S.pathAxes(ip), S.convAxes(ip), S.last, alg);

            if drawn(ip)
                ylims(end+1, :) = yl; %#ok<AGROW>
            else
                cla(S.pathAxes(ip), 'reset');  cla(S.convAxes(ip), 'reset');
                axis(S.pathAxes(ip), 'off');   axis(S.convAxes(ip), 'off');
                title(S.pathAxes(ip), sprintf( ...
                    ['No %s refinement history in this run.\n\n' ...
                     'Tick "record refinement steps", enable %s, and keep the\n' ...
                     'traced atom index at or below the scatterer count.'], ...
                    upper(alg), upper(alg)), ...
                    'FontWeight', 'normal', 'FontSize', 11);
            end
        end

        % a shared y-range makes the two convergence curves comparable at a
        % glance, which is the whole point of stacking them
        if size(ylims, 1) == 2
            shared = [min(ylims(:,1)), max(ylims(:,2))];
            for ip = 1:numel(S.pathAlgs)
                if drawn(ip)
                    ylim(S.convAxes(ip), shared);
                end
            end
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

        out = S.last;

        % pop out what the user is looking at: the images tab gives the
        % images, the positions tab the positions, anything else gives both
        selected = S.tabs.SelectedTab;
        if isequal(selected, S.imageTab)
            panels = result_panels(out.x_hat, 'image');
        elseif isequal(selected, S.positionTab)
            panels = result_panels(out.x_hat, 'positions');
        else
            panels = result_panels(out.x_hat);
        end
        n = size(panels, 1);

        if n == 0
            return
        end

        [rows, cols] = panel_geometry(n);
        f = figure('Name', scenario_title(out.cfg), ...
            'Position', [100 100 520*cols 520*rows]);
        for i = 1:n
            ax = subplot(rows, cols, i, 'Parent', f);
            draw_panel(ax, panels{i,1}, panels{i,2}, out, ...
                S.logScaleCheck.Value, S.dynRangeField.Value, ...
                S.ambViewCheck.Value, S.aliasCheck.Value);
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

            summary(i, :) = { alg_label(alg), r.error, size(r.pairs, 1), ...
                numel(r.missed), numel(r.false_alarms), ...
                median(r.d), max(r.d) };

            for j = 1:size(r.pairs, 1)
                it = r.pairs(j, 1);   % index into the true locations
                ie = r.pairs(j, 2);   % index into the estimated locations
                pairs(end+1, :) = { alg_label(alg), ...
                    out.target_locations(it, 1), out.target_locations(it, 2), ...
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
            err('omp'), err('mod_omp'), err('nomp'), err('promp'), err('bp'), ...
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
    for a = {'omp', 'mod_omp', 'nomp', 'promp', 'bp'}
        if is_scored(x_hat, a{1})
            algs{end+1} = a{1}; %#ok<AGROW>
        end
    end
end

function label = alg_label(alg)
% ALG_LABEL  Display name for an algorithm key.
%   Axes titles interpret '_' as a TeX subscript, so 'mod_omp' would render
%   as MOD with a subscripted O.

    switch alg
        case 'mod_omp'
            label = 'Mod-OMP';
        otherwise
            label = upper(alg);
    end
end

function tf = is_scored(x_hat, alg)
    tf = isfield(x_hat, alg) && isfield(x_hat.(alg), 'error');
end

function panels = result_panels(x_hat, want_mode)
% RESULT_PANELS  The plots to show, as {algorithm, mode} rows.
%
%   An algorithm can contribute more than one view. OMP gets both: the
%   gridded image, and the extracted peak positions drawn the same way as
%   NOMP's, which is the easier of the two to compare against truth when
%   the scatterers are closely spaced.
%
%   RESULT_PANELS(x_hat, mode) keeps only the rows of that mode, which is how
%   the images and the positions end up in separate tabs.

    if nargin < 2
        want_mode = '';
    end

    panels = cell(0, 2);

    if is_scored(x_hat, 'omp')
        panels(end+1, :) = {'omp', 'image'};
        panels(end+1, :) = {'omp', 'positions'};
    end

    if is_scored(x_hat, 'mod_omp')
        panels(end+1, :) = {'mod_omp', 'image'};
        panels(end+1, :) = {'mod_omp', 'positions'};
    end

    if is_scored(x_hat, 'nomp')
        panels(end+1, :) = {'nomp', 'positions'};
    end

    if is_scored(x_hat, 'promp')
        panels(end+1, :) = {'promp', 'positions'};
    end

    if is_scored(x_hat, 'bp')
        panels(end+1, :) = {'bp', 'image'};
    end

    if ~isempty(want_mode) && ~isempty(panels)
        panels = panels(strcmp(panels(:,2), want_mode), :);
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
    cols = {'Run', 'RMS_OMP_m', 'RMS_ModOMP_m', 'RMS_NOMP_m', 'RMS_PROMP_m', 'RMS_BP_m', ...
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

function [ok, ylims] = draw_refinement_path(ax_trace, ax_conv, out, alg)
% DRAW_REFINEMENT_PATH  One algorithm's per-step position estimates for a
% single atom, as a chain of arrows zoomed on the true scatterer it converged
% to.
%
%   ALG is 'nomp' or 'promp'. The two record different things and the plot is
%   the same either way: NOMP's history is its per-atom Newton steps plus each
%   later cyclic refinement of that atom, PROMP's is the atom's position after
%   every joint Gauss-Newton step it takes part in, which continues on every
%   OMP iteration after the one that selected it.
%
%   Left axes  -- the walk itself: a square marks the grid node OMP selected,
%                 arrows step from each estimate to the next, coloured light
%                 to dark with iteration, and a red circle is the truth. The
%                 view is centred on the truth so the residual offset is what
%                 you read off the axes, not the absolute position.
%   Right axes -- distance to that truth per step, on a log scale, with the
%                 pixel pitch drawn in so you can see when the estimate goes
%                 sub-pixel.
%
%   Returns false when the run carries no history, so the caller can put up
%   an explanation instead of an empty pair of axes.

    ok    = false;
    ylims = [0 1];

    if ~isfield(out.x_hat, alg) || ~isfield(out.x_hat.(alg), 'p_hat_hist')
        return
    end

    H = out.x_hat.(alg).p_hat_hist;         % [nsteps x 2], columns [x, y]
    if isempty(H) || size(H, 2) ~= 2
        return
    end

    % the history is trimmed at the source to the steps actually taken, so
    % every row here is a real estimate -- including a legitimate [0 0] if the
    % detection step happened to pick the centre pixel
    if size(H, 1) < 2
        return
    end

    u0 = out.u0;

    % which atom was traced, and the truth it ended up closest to
    idx = 1;
    if isfield(out.opts, 'save_atom_idx')
        idx = out.opts.save_atom_idx;
    end
    idx = max(1, min(idx, size(out.x_hat.(alg).positions, 1)));
    final_est = out.x_hat.(alg).positions(idx, :);

    T = out.latent_locations;
    [~, it] = min((T(:,1) - final_est(1)).^2 + (T(:,2) - final_est(2)).^2);
    tx = T(it, 1);
    ty = T(it, 2);

    n    = size(H, 1);
    cmap = parula(max(n - 1, 1));

    % ---------------------------------------------------- the walk
    cla(ax_trace, 'reset');
    hold(ax_trace, 'on');

    for i = 1:n-1
        quiver(ax_trace, H(i,1), H(i,2) + u0, ...
            H(i+1,1) - H(i,1), H(i+1,2) - H(i,2), 0, ...
            'Color', cmap(i,:), 'LineWidth', 1.6, 'MaxHeadSize', 0.5);
    end

    plot(ax_trace, H(:,1), H(:,2) + u0, 'o', 'MarkerSize', 4, ...
        'MarkerFaceColor', 'w', 'MarkerEdgeColor', [0.4 0.4 0.4]);
    plot(ax_trace, H(1,1), H(1,2) + u0, 's', 'MarkerSize', 13, ...
        'LineWidth', 1.6, 'MarkerEdgeColor', [0 0.45 0.74]);
    plot(ax_trace, H(end,1), H(end,2) + u0, 'p', 'MarkerSize', 15, ...
        'LineWidth', 1.2, 'MarkerFaceColor', [0.20 0.65 0.35], ...
        'MarkerEdgeColor', 'k');
    plot(ax_trace, tx, ty + u0, 'o', 'MarkerSize', 13, ...
        'LineWidth', 1.8, 'MarkerEdgeColor', [1 0 0]);

    % label every step when there are few, otherwise thin them out -- a long
    % trace (one entry per cyclic refinement of this atom, so it grows with
    % the scatterer count) turns into unreadable overprinted text otherwise
    if n <= 12
        label_at = 1:n;
    else
        label_at = unique([1, 1:ceil(n/8):n, n]);
    end
    for i = label_at
        text(ax_trace, H(i,1), H(i,2) + u0, sprintf('  %d', i - 1), ...
            'FontSize', 8, 'Color', [0.35 0.35 0.35], 'Parent', ax_trace);
    end

    % centre on the truth, with a floor so a converged trace does not zoom
    % into numerical noise
    span  = max([abs(H(:,1) - tx); abs(H(:,2) - ty)]);
    pitch = max(out.sim_config.CrossRangePixelRes_m, ...
                out.sim_config.RangePixelRes_m);
    half  = max(span * 1.35, 0.02 * pitch);

    xlim(ax_trace, tx + [-half, half]);
    ylim(ax_trace, ty + u0 + [-half, half]);
    set(ax_trace, 'YDir', 'reverse');
    grid(ax_trace, 'on');
    box(ax_trace, 'on');
    xlabel(ax_trace, 'crossrange [m]');
    ylabel(ax_trace, 'range [m]');
    title(ax_trace, sprintf('%s atom %d: %d steps, final miss %.3g m', ...
        upper(alg), idx, n - 1, hypot(H(end,1) - tx, H(end,2) - ty)));
    hold(ax_trace, 'off');

    % ---------------------------------------------- convergence curve
    d = hypot(H(:,1) - tx, H(:,2) - ty);

    % a step can land exactly on the truth, and log(0) would stretch the axis
    % over hundreds of decades. sub-nanometre is numerically converged for a
    % 0.3 m wavelength, so that is the floor and the axis is pinned to it
    floor_d = 1e-9;
    d_plot  = max(d, floor_d);

    cla(ax_conv, 'reset');
    semilogy(ax_conv, 0:n-1, d_plot, '-o', ...
        'Color', [0 0.45 0.74], 'LineWidth', 1.4, 'MarkerSize', 4, ...
        'MarkerFaceColor', [0 0.45 0.74]);
    hold(ax_conv, 'on');
    yline(ax_conv, out.sim_config.CrossRangePixelRes_m, '--', ...
        'crossrange pixel', 'Color', [0.6 0.6 0.6], 'FontSize', 8);
    yline(ax_conv, out.sim_config.RangePixelRes_m, '--', ...
        'range pixel', 'Color', [0.6 0.6 0.6], 'FontSize', 8);
    grid(ax_conv, 'on');
    box(ax_conv, 'on');
    xlabel(ax_conv, 'refinement step');
    ylabel(ax_conv, 'distance to truth [m]');
    title(ax_conv, sprintf('%s convergence', upper(alg)));
    xlim(ax_conv, [-0.5, n - 0.5]);
    ylim(ax_conv, [floor_d / 2, max(d_plot) * 2]);
    hold(ax_conv, 'off');

    ylims = ylim(ax_conv);
    ok    = true;
end

function [v, base, vlabel] = one_dim_values(v, log_scale, dyn_range)
% ONE_DIM_VALUES  Amplitudes for a crossrange stem plot, linear or in dB.
%
%   In dB the trace is floored at DYN_RANGE below the peak, so the stems hang
%   from a sensible baseline instead of running off to -Inf wherever the
%   reconstruction is exactly zero -- which, for OMP, is almost everywhere.

    if log_scale
        v   = 20*log10(max(v, realmin));
        top = max(v(isfinite(v)));
        if isempty(top)
            top = 0;
        end
        base   = top - dyn_range;
        v      = max(v, base);
        vlabel = 'Log-Scaled [dB]';
    else
        base   = 0;
        vlabel = 'Amplitude [Linear]';
    end
end

function draw_ambiguity_bands(ax, Wx)
% DRAW_AMBIGUITY_BANDS  Mark where one crossrange ambiguity ends and the next
% begins, and name each band.
%
%   Crossrange repeats with period Wx, so band n spans [(n-0.5)*Wx,
%   (n+0.5)*Wx]. The dividers are only drawn when the axes actually show more
%   than one band -- on a single-band view they would just be two lines at the
%   edges of the plot.

    if isempty(Wx) || ~isfinite(Wx) || Wx <= 0
        return
    end

    xl = xlim(ax);
    yl = ylim(ax);

    if diff(xl) < 1.2 * Wx
        return
    end

    wasHold = ishold(ax);
    hold(ax, 'on');

    for n = -3:3
        edge = (n + 0.5) * Wx;
        if edge > xl(1) && edge < xl(2)
            plot(ax, [edge edge], yl, '--', 'Color', [0.55 0.55 0.55], ...
                'LineWidth', 1, 'HandleVisibility', 'off');
        end

        centre = n * Wx;
        if centre > xl(1) && centre < xl(2)
            text(ax, centre, yl(1), sprintf(' amb %+d ', n), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'top', ...
                'FontSize', 8, 'FontWeight', 'bold', ...
                'Color', [0.35 0.35 0.35], 'Parent', ax);
        end
    end

    if ~wasHold
        hold(ax, 'off');
    end
end

function draw_panel(ax, alg, mode, out, log_scale, dyn_range, show_ambiguities, show_aliased)
% DRAW_PANEL  Render one view of one algorithm's result into AX.
%
%   mode 'image'     -- the gridded latent image, truth overlaid
%   mode 'positions' -- the recovered scatterer positions alone, truth
%                       overlaid, identically for OMP and NOMP
%
%   SHOW_AMBIGUITIES widens the view to the three ambiguity bands -1, 0 and +1
%   and, for image panels, repeats the formed image at -Wx and +Wx. The image
%   former only ever solves for the bands it was given, but crossrange
%   genuinely repeats with period Wx, so those replicas are the other
%   positions that explain the same measurement equally well. Off, the panel
%   shows only what was formed.
%
%   SHOW_ALIASED draws the magenta markers at the folded positions of the
%   ambiguous scatterers. Useful for reading ghosts inside the formed band,
%   and worth switching off once the point is made -- or when SHOW_AMBIGUITIES
%   is on and the true positions are already visible in their own bands.
%
%   Works with both uiaxes (in the GUI) and ordinary axes (in a popped out
%   figure), so the two always show the same thing.

    if nargin < 7 || isempty(show_ambiguities)
        show_ambiguities = false;
    end
    if nargin < 8 || isempty(show_aliased)
        show_aliased = true;
    end

    Wx = [];
    if isfield(out, 'sim_config') && isfield(out.sim_config, 'W_x_m')
        Wx = out.sim_config.W_x_m;
    end
    show_ambiguities = show_ambiguities && ~isempty(Wx);

    y_array = out.y_array;
    u0      = out.u0;
    r       = out.x_hat.(alg);

    % The modified OMP forms one image per ambiguity, stacked along
    % crossrange, so its panel carries a wider axis of its own and already
    % covers the bands the other algorithms can only show as replicas.
    stacked = isfield(r, 'x_array') && ~isempty(r.x_array);
    if stacked
        x_array = r.x_array;
    else
        x_array = out.x_array;
    end

    % with the wider view the neighbouring bands are on screen, so every
    % scatterer can be drawn where it actually is rather than only the ones
    % the image former covers
    if show_ambiguities || stacked
        truth = out.target_locations;
    else
        truth = out.latent_locations;
    end
    true_x = truth(:, 1);
    true_y = truth(:, 2) + u0;

    % Crossrange is only unambiguous over a width Wx. A scatterer sitting in
    % ambiguity n has its energy folded to x - n*Wx, so it shows up in the
    % image at that folded position with nothing underneath it. Drawing those
    % folded positions is what tells a ghost apart from a real scatterer --
    % especially when the image former spans fewer ambiguities than the target
    % occupies, in which case the true position is off the plot entirely.
    alias_x = [];
    alias_y = [];
    if show_aliased && isfield(out, 'amb_of_k') && isfield(out.sim_config, 'W_x_m')
        Wx  = out.sim_config.W_x_m;
        amb = out.amb_of_k(:);
        n   = min(numel(amb), size(out.target_locations, 1));
        sel = false(n, 1);
        sel(1:n) = amb(1:n) ~= 0;      % ambiguity 0 folds onto itself

        if any(sel)
            alias_x = out.target_locations(sel, 1) - amb(sel) * Wx;
            alias_y = out.target_locations(sel, 2) + u0;
        end
    end

    % one range cell means y_array is a single value: there is no range axis
    % to plot against, so every panel becomes a crossrange line
    is_range_cell = isscalar(y_array);

    cla(ax);

    switch mode
        case 'image'
            img = abs(r.image);

            if is_range_cell
                [v, base, vlabel] = one_dim_values(img(:).', log_scale, dyn_range);
                h = stem(ax, x_array, v, 'filled', 'MarkerSize', 3, ...
                    'Color', [0 0.45 0.74], 'LineWidth', 1);
                h.BaseValue = base;
                ylim(ax, [base, max([v, base + eps]) * 1.05]);
                ylabel(ax, vlabel);
                heading = alg_label(alg);

            elseif log_scale
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

            if ~is_range_cell

            % the formed image covers its own bands only; draw it again either
            % side so the ambiguous replicas sit next to the true positions.
            % A stacked image already spans them, so it needs no replicas.
            if show_ambiguities && ~stacked
                hold(ax, 'on');
                for off = [-Wx, Wx]
                    imagesc(ax, x_array + off, y_array + u0, img);
                end
                hold(ax, 'off');
            end

                colormap(ax, gray);
                cb = colorbar(ax);
                cb.Label.String = cbar_label;
                heading = alg_label(alg);
            end

        case 'positions'
            % OMP's positions are the extracted image peaks and NOMP's are
            % estimated off-grid directly, but both are [x, y] in metres, so
            % they are drawn the same way and read the same way
            if is_range_cell
                % every estimate is in the one cell, so the informative axis
                % is the recovered amplitude against crossrange
                amp = ones(size(r.positions, 1), 1);
                if isfield(r, 'alpha') && numel(r.alpha) == size(r.positions, 1)
                    amp = abs(r.alpha(:));
                end
                stem(ax, r.positions(:, 1), amp, 'filled', 'MarkerSize', 5, ...
                    'Color', [0 0 1], 'LineWidth', 1.2);
                ylabel(ax, 'Reflectivity magnitude');
            else
                plot(ax, r.positions(:, 1), r.positions(:, 2) + u0, '*', ...
                    'MarkerEdgeColor', [0 0 1], 'MarkerSize', 10, ...
                    'LineWidth', 1.5);
                set(ax, 'YDir', 'reverse');
            end
            heading = [alg_label(alg) ' positions'];

        otherwise
            error('isar_gui:panelMode', 'unknown panel mode ''%s''', mode);
    end

    hold(ax, 'on');
    if is_range_cell
        % with amplitude on the vertical axis a marker has nowhere sensible to
        % sit, so the truth becomes a line at its crossrange
        h_truth = xline(ax, true_x(1), '-', 'Color', [1 0 0], 'LineWidth', 1.5);
        for it = 2:numel(true_x)
            xline(ax, true_x(it), '-', 'Color', [1 0 0], 'LineWidth', 1.5, ...
                'HandleVisibility', 'off');
        end
    else
        h_truth = plot(ax, true_x, true_y, 'o', 'MarkerEdgeColor', [1 0 0], ...
            'MarkerSize', 10, 'LineWidth', 1.5);
    end

    if ~isempty(alias_x)
        if is_range_cell
            h_alias = xline(ax, alias_x(1), '-', 'Color', [1 0 1], 'LineWidth', 1.5);
            for ia = 2:numel(alias_x)
                xline(ax, alias_x(ia), '-', 'Color', [1 0 1], 'LineWidth', 1.5, ...
                    'HandleVisibility', 'off');
            end
        else
            h_alias = plot(ax, alias_x, alias_y, 'd', 'MarkerEdgeColor', [1 0 1], ...
                'MarkerSize', 10, 'LineWidth', 1.5);
        end
        legend(ax, [h_truth, h_alias], ...
            {'true scatterer', 'aliased into ambiguity 0'}, ...
            'Location', 'southoutside', 'Orientation', 'horizontal', ...
            'FontSize', 8, 'Box', 'off');
    end
    hold(ax, 'off');

    if ~is_range_cell
        axis(ax, 'square');
    end
    if show_ambiguities && ~stacked
        xlim(ax, [-1.5, 1.5] * Wx);
    else
        xlim(ax, [min(x_array), max(x_array)]);
    end
    if ~is_range_cell
        ylim(ax, [min(y_array) + u0, max(y_array) + u0]);
    end

    draw_ambiguity_bands(ax, Wx);

    xlabel(ax, 'crossrange [m]');
    if ~is_range_cell
        ylabel(ax, 'range [m]');
    end
    title(ax, sprintf('%s  (RMS error: %.3f m)', heading, r.error));
end
