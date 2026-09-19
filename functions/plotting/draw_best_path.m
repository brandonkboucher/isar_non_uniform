function draw_best_path(ax, pp, n_acc, color, u0, r, ca, fs, min_label_dist)
% DRAW_BEST_PATH  Highlight one Newton trajectory: accepted iterates as a
% thick outlined line with markers, a rejected final step dashed and ending
% in an x. With MIN_LABEL_DIST > 0 each iterate is numbered with its
% correlation |a'r|/(||a|| ||r||); iterates closer than MIN_LABEL_DIST to the
% previous labelled one are skipped so converged steps do not pile up.

    y = pp(2,:) + u0;
    plot(ax, pp(1,1:n_acc), y(1:n_acc), '-', 'Color', 'k', 'LineWidth', 4.5, ...
        'HandleVisibility', 'off');
    plot(ax, pp(1,1:n_acc), y(1:n_acc), '-o', 'Color', color, 'LineWidth', 2.5, ...
        'MarkerSize', 6, 'MarkerFaceColor', color, 'MarkerEdgeColor', 'k', ...
        'HandleVisibility', 'off');
    if n_acc < size(pp,2)
        plot(ax, pp(1,n_acc:end), y(n_acc:end), '--', 'Color', color, 'LineWidth', 2, ...
            'HandleVisibility', 'off');
        plot(ax, pp(1,end), y(end), 'x', 'MarkerSize', 12, 'LineWidth', 2.5, ...
            'Color', color, 'HandleVisibility', 'off');
    end

    if min_label_dist <= 0
        return
    end
    last = [Inf; Inf];
    for s = 1:size(pp,2)
        is_last_acc = s == n_acc;
        is_rej      = s > n_acc;
        if norm(pp(:,s) - last) < min_label_dist && ~is_last_acc && ~is_rej
            continue
        end
        a  = ca(pp(1,s), pp(2,s));
        cr = abs(a' * r) / (norm(a) * norm(r));
        lbl = sprintf(' %d: %.3f', s - 1, cr);
        if is_rej, lbl = [lbl ' (rejected)']; end %#ok<AGROW>
        text(ax, pp(1,s), y(s), lbl, 'Color', color, 'FontSize', fs - 2, ...
            'FontWeight', 'bold', 'BackgroundColor', 'k', 'Margin', 1, ...
            'VerticalAlignment', 'bottom', 'Clipping', 'on');
        last = pp(:,s);
    end
end
