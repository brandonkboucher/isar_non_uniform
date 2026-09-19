function plot_mod_omp_candidates(cand, truth, grd, u0, theta_m, f_hat_l, fc, sc, options)
% PLOT_MOD_OMP_CANDIDATES  One figure per mod-OMP iteration showing where the
% doppelganger candidates start, where Newton takes them, and which one wins.
%
%   Top: |a(x,y)' r| / (||a|| ||r||) over every ambiguity band for that
%   iteration's residual, with the truth, the band-0 detection, the seeds and
%   the Newton end points. Bottom: a fine map around each ambiguity's search
%   window with the Newton paths, and the end-point correlation of every
%   candidate against its seed offset, which is the contest mod-OMP decides.
%
%   An iteration is attributed to the true scatterer whose folded position
%   (crossrange reduced into band 0) is nearest the detection, since what the
%   one-ambiguity image former detects is that scatterer's ghost.
%
%   The correlation maps are built with BUILD_AS, one atom per pixel.
%
%   options.plot_best_candidate_paths highlights, for every ambiguity, the
%   full Newton trajectory of that ambiguity's best candidate: seed, each
%   iterate (numbered, with its correlation in the zoom panels) and any
%   rejected final step.
%
%   See also BEST_CANDIDATE_PATH, DRAW_BEST_PATH, BUILD_AS, MOD_OMP_NEWTON.

    ura       = options.use_range_approx;
    ca        = @(x, y) compute_atom(x, y, u0, theta_m, f_hat_l, fc, ura);
    show_best = isfield(options, 'plot_best_candidate_paths') ...
        && options.plot_best_candidate_paths;
    Wx        = grd.Wx;
    res_x     = grd.cross_range_pixel_res;
    res_y     = grd.range_pixel_res;
    amb_index = unique(cand(1).amb);
    n_amb     = numel(amb_index);
    colors    = lines(n_amb);
    fs        = 12;

    band   = round(truth(:,1) / Wx);
    folded = [truth(:,1) - band*Wx, truth(:,2)];

    plot_all = isfield(options, 'plot_candidates_all_iterations') ...
        && options.plot_candidates_all_iterations;

    % which iterations to draw
    k_of_iter = zeros(numel(cand), 1);
    for i = 1:numel(cand)
        [~, k_of_iter(i)] = min(vecnorm(folded - [cand(i).x0, cand(i).y0], 2, 2));
    end
    iters = find(plot_all | band(k_of_iter) ~= 0).';
    if isempty(iters)
        fprintf('plot_mod_omp_candidates: no iteration detected an aliased scatterer\n');
        return
    end

    % overview dictionary spanning every searched band, built once
    x_ov = (min(amb_index) - 0.5)*Wx : res_x : (max(amb_index) + 0.5)*Wx;
    y_ov = grd.y_array;
    [Xo, Yo] = meshgrid(x_ov, y_ov);
    A_ov = build_AS([Xo(:).'; Yo(:).'], u0, theta_m, f_hat_l, fc, ura);

    for i = iters

        c  = cand(i);
        r  = c.r;
        k  = k_of_iter(i);
        ib = c.best;
        pb = c.p(:, ib);

        n_rows = 2 + show_best;
        f = figure('Visible', 'off', 'Position', [0 0 1900 525*n_rows]);
        tl = tiledlayout(f, n_rows, n_amb + 1, 'TileSpacing', 'compact', 'Padding', 'compact');

        %------------------------- overview -----------------------------
        ax = nexttile(tl, [1, n_amb + 1]);
        map = reshape(abs(A_ov' * r) ./ (vecnorm(A_ov).' * norm(r)), size(Xo));
        imagesc(ax, x_ov, y_ov + u0, map); colormap(ax, gray); clim(ax, [0 1]);
        cb = colorbar(ax); cb.Label.String = '|a^H r| / (||a|| ||r||)';
        hold(ax, 'on')
        for edge = ((min(amb_index) - 0.5):(max(amb_index) + 0.5)) * Wx
            plot(ax, [edge edge], [y_ov(1) y_ov(end)] + u0, '--', ...
                'Color', [1 1 1]*0.8, 'HandleVisibility', 'off');
        end
        h = gobjects(0);
        h(end+1) = plot(ax, truth(:,1), truth(:,2) + u0, 'o', 'MarkerSize', 12, ...
            'LineWidth', 2, 'Color', [0.9 0.1 0.1], 'DisplayName', 'truth');  %#ok<AGROW>
        for kk = 1:size(truth,1)
            text(ax, truth(kk,1), truth(kk,2) + u0, sprintf('  %d (band %+d)', kk, band(kk)), ...
                'Color', [0.9 0.1 0.1], 'FontSize', fs, 'FontWeight', 'bold');
        end
        h(end+1) = plot(ax, c.x0, c.y0 + u0, 'd', 'MarkerSize', 12, 'LineWidth', 2, ...
            'Color', [1 0.8 0], 'DisplayName', 'detection (band-0 grid)');  %#ok<AGROW>
        for ia = 1:n_amb
            m = c.amb == amb_index(ia);
            h(end+1) = plot(ax, c.p(1,m), c.p(2,m) + u0, 'x', 'MarkerSize', 8, 'LineWidth', 1.5, ...
                'Color', colors(ia,:), 'DisplayName', sprintf('Newton ends, amb %+d', amb_index(ia))); %#ok<AGROW>
        end
        h(end+1) = plot(ax, pb(1), pb(2) + u0, 'p', 'MarkerSize', 20, 'LineWidth', 2, ...
            'Color', [0 1 0.4], 'DisplayName', 'selected');  %#ok<AGROW>
        hold(ax, 'off')
        set(ax, 'YDir', 'reverse', 'FontSize', fs)
        xlabel(ax, 'crossrange [m]'); ylabel(ax, 'range [m]')
        legend(ax, h, 'Location', 'eastoutside', 'FontSize', fs - 1)
        title(ax, sprintf('residual correlation at iteration %d, over %d ambiguities', i, n_amb))

        %---------------------- zoom per ambiguity ----------------------
        for ia = 1:n_amb
            ax = nexttile(tl);
            m  = c.amb == amb_index(ia);
            xc = c.x0 + amb_index(ia)*Wx;
            xl = xc + [-1 1] * (sc.num_offsets_pixels + 4) * res_x;

            pts = [c.seed(:,m), c.p(:,m), cell2mat(c.path(m))];
            yl  = [min([pts(2,:), c.y0 - 0.75]) - 0.2, max([pts(2,:), c.y0 + 0.75]) + 0.2];

            xz = xl(1):res_x/3:xl(2);
            yz = yl(1):res_y/3:yl(2);
            [Xz, Yz] = meshgrid(xz, yz);
            Az = build_AS([Xz(:).'; Yz(:).'], u0, theta_m, f_hat_l, fc, ura);
            mz = reshape(abs(Az' * r) ./ (vecnorm(Az).' * norm(r)), size(Xz));

            imagesc(ax, xz, yz + u0, mz); colormap(ax, gray); clim(ax, [0 1]);
            hold(ax, 'on')
            % Newton's history records a step before the acceptance test, so
            % a final iterate that differs from where the search ended was
            % rejected: draw the accepted path solid and that step dotted
            for ic = find(m)
                pp = [c.seed(:,ic), c.path{ic}];
                rejected = size(pp,2) > 1 && any(pp(:,end) ~= c.p(:,ic));
                n_acc = size(pp,2) - rejected;
                plot(ax, pp(1,1:n_acc), pp(2,1:n_acc) + u0, '-', 'Color', colors(ia,:), 'LineWidth', 1.2);
                if rejected
                    plot(ax, pp(1,n_acc:end), pp(2,n_acc:end) + u0, ':', 'Color', colors(ia,:), 'LineWidth', 1.2);
                    plot(ax, pp(1,end), pp(2,end) + u0, 'o', 'MarkerSize', 5, 'Color', colors(ia,:));
                end
            end
            plot(ax, c.seed(1,m), c.seed(2,m) + u0, '.', 'MarkerSize', 12, 'Color', colors(ia,:));
            plot(ax, c.p(1,m), c.p(2,m) + u0, 'x', 'MarkerSize', 8, 'LineWidth', 1.5, 'Color', colors(ia,:));
            in = truth(:,1) >= xl(1) & truth(:,1) <= xl(2);
            plot(ax, truth(in,1), truth(in,2) + u0, 'o', 'MarkerSize', 14, 'LineWidth', 2, 'Color', [0.9 0.1 0.1]);
            if c.amb(ib) == amb_index(ia)
                plot(ax, pb(1), pb(2) + u0, 'p', 'MarkerSize', 20, 'LineWidth', 2, 'Color', [0 1 0.4]);
            end
            hold(ax, 'off')
            set(ax, 'YDir', 'reverse', 'FontSize', fs)
            xlim(ax, xl); ylim(ax, yl + u0)
            xlabel(ax, 'crossrange [m]'); ylabel(ax, 'range [m]')
            stuck = sum(all(c.p(:,m) == c.seed(:,m), 1));
            title(ax, sprintf('ambiguity %+d: best %.3f, %d/%d never left seed', ...
                amb_index(ia), max(c.corr(m)), stuck, nnz(m)))
        end

        %------------------- the contest itself -------------------------
        ax = nexttile(tl);
        hold(ax, 'on')
        for ia = 1:n_amb
            m   = c.amb == amb_index(ia);
            off = c.seed(1,m) - (c.x0 + amb_index(ia)*Wx);
            plot(ax, off, c.corr(m), '-o', 'Color', colors(ia,:), 'LineWidth', 1.5, ...
                'MarkerSize', 4, 'DisplayName', sprintf('amb %+d', amb_index(ia)));
        end
        plot(ax, c.seed(1,ib) - (c.x0 + c.amb(ib)*Wx), c.corr(ib), 'p', 'MarkerSize', 18, ...
            'LineWidth', 2, 'Color', [0 0.6 0.25], 'DisplayName', 'selected');
        hold(ax, 'off'); grid(ax, 'on'); box(ax, 'on')
        set(ax, 'FontSize', fs)
        xlabel(ax, 'seed offset from x_0 + k W_x [m]')
        ylabel(ax, 'correlation after Newton')
        legend(ax, 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', fs - 1)
        title(ax, 'candidate contest')

        %-------------- best candidate's trajectory, per ambiguity --------
        % A tight window around each ambiguity's best candidate, so steps of a
        % few millimetres are visible, then its correlation against step.
        if show_best
            best_corr = cell(1, n_amb);
            best_nacc = zeros(1, n_amb);
            for ia = 1:n_amb
                ax = nexttile(tl);
                [pp, n_acc] = best_candidate_path(c, amb_index(ia));

                % window: the path's bounding box, padded, never smaller than
                % a quarter pixel either side
                half = max((max(pp,[],2) - min(pp,[],2)) / 2 * 1.4, [res_x; res_y] / 4);
                ctr  = (max(pp,[],2) + min(pp,[],2)) / 2;
                xw = ctr(1) + [-1 1]*half(1);
                yw = ctr(2) + [-1 1]*half(2);

                xz = linspace(xw(1), xw(2), 61);
                yz = linspace(yw(1), yw(2), 61);
                [Xz, Yz] = meshgrid(xz, yz);
                Az = build_AS([Xz(:).'; Yz(:).'], u0, theta_m, f_hat_l, fc, ura);
                mz = reshape(abs(Az' * r) ./ (vecnorm(Az).' * norm(r)), size(Xz));

                imagesc(ax, xz, yz + u0, mz); colormap(ax, gray); clim(ax, [0 1]);
                hold(ax, 'on')
                in = truth(:,1) >= xw(1) & truth(:,1) <= xw(2) ...
                   & truth(:,2) >= yw(1) & truth(:,2) <= yw(2);
                plot(ax, truth(in,1), truth(in,2) + u0, 'o', 'MarkerSize', 14, ...
                    'LineWidth', 2, 'Color', [0.9 0.1 0.1]);
                draw_best_path(ax, pp, n_acc, colors(ia,:), u0, r, ca, fs, ...
                    0.04 * max(diff(xw), diff(yw)));
                hold(ax, 'off')
                set(ax, 'YDir', 'reverse', 'FontSize', fs)
                xlim(ax, xw); ylim(ax, yw + u0)
                xlabel(ax, 'crossrange [m]'); ylabel(ax, 'range [m]')

                % correlation at each point of the path, for the last tile
                cr = zeros(1, size(pp,2));
                for s = 1:size(pp,2)
                    a = ca(pp(1,s), pp(2,s));
                    cr(s) = abs(a' * r) / (norm(a) * norm(r));
                end
                best_corr{ia} = cr;
                best_nacc(ia) = n_acc;

                title(ax, sprintf('amb %+d best candidate: %d accepted step(s)%s\n%.3f at seed -> %.3f', ...
                    amb_index(ia), n_acc - 1, repmat(', 1 rejected', 1, n_acc < size(pp,2)), ...
                    cr(1), cr(n_acc)))
            end

            % correlation along each best trajectory; step 0 is the seed
            ax = nexttile(tl);
            hold(ax, 'on')
            for ia = 1:n_amb
                cr = best_corr{ia}; n_acc = best_nacc(ia);
                plot(ax, 0:n_acc-1, cr(1:n_acc), '-o', 'Color', colors(ia,:), 'LineWidth', 1.8, ...
                    'MarkerFaceColor', colors(ia,:), 'DisplayName', sprintf('amb %+d', amb_index(ia)));
                if n_acc < numel(cr)
                    plot(ax, [n_acc-1 n_acc], cr(n_acc:n_acc+1), '--o', 'Color', colors(ia,:), ...
                        'LineWidth', 1.5, 'MarkerFaceColor', 'w', 'HandleVisibility', 'off');
                end
            end
            hold(ax, 'off'); grid(ax, 'on'); box(ax, 'on')
            set(ax, 'FontSize', fs)
            xlabel(ax, 'Newton step (0 = seed)')
            ylabel(ax, 'correlation')
            legend(ax, 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', fs - 1)
            title(ax, 'best candidate per ambiguity (dashed, hollow = rejected step)')
        end

        sgtitle(tl, sprintf(['iteration %d: detected (%.2f, %.2f) on the band-0 grid  |  ' ...
            'ghost of truth %d at (%.2f, %.2f), band %+d  |  selected band %+d at (%.2f, %.2f)'], ...
            i, c.x0, c.y0, k, truth(k,1), truth(k,2), band(k), c.amb(ib), pb(1), pb(2)), ...
            'FontSize', fs + 2, 'FontWeight', 'bold');

        if options.save_plots
            saveas(f, fullfile('plots', sprintf('mod_omp_candidates_iter%02d.png', i)));
            close(f)
        else
            f.Visible = 'on';
        end
    end
end
