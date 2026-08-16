function [outputArg1,outputArg2] = plot_all_pairings(...
    scenario_number, ...
    target_locations,...
    u0,...
    x_array, ...
    y_array, ...
    x_hat, ...
    options, ...
    super_title)

    num_plots = 0; 
    font_size = 24;

    if isfield(options, 'execute_itsa') ...
            && options.execute_itsa
        num_plots = num_plots + 1;
    end

    if isfield(options, 'execute_omp') ...
            && options.execute_omp
        num_plots = num_plots + 1;
    end

    if isfield(options, 'execute_nomp') ...
            && options.execute_nomp
        num_plots = num_plots + 1;
    end

    if isfield(options, 'execute_sbl') ...
            && options.execute_sbl
        num_plots = num_plots + 1;
    end

    if isfield(options, 'execute_bp') ...
            && options.execute_bp
        num_plots = num_plots + 1;
    end

    f = figure('Visible','off');
    iplot = 1;

    if isfield(options, 'execute_nomp') ...
            && options.execute_nomp

        T = target_locations; 
        E = x_hat.nomp.positions;
        d = x_hat.nomp.d;
        pairs = x_hat.nomp.pairs;

        subplot(1,num_plots,iplot)
        iplot = iplot + 1;

        T(:,2) = T(:,2) + u0;
        E(:,2) = E(:,2) + u0;
    
        if ~isempty(pairs)
            if max(pairs(:,1)) > size(T,1) || max(pairs(:,2)) > size(E,1) ...
                    || min(pairs(:)) < 1
                error('pairs contains indices outside the location arrays.');
            end
        end
        
        h = gca;
    
        uT = setdiff((1:size(T,1)).', pairs(:,1));   % unmatched true  -> missed
        uE = setdiff((1:size(E,1)).', pairs(:,2));   % unmatched est.  -> false alarm
    
        wasHold = ishold(h);
        hold(h, 'on');
    
        % connecting segments first so markers draw on top
        for k = 1:size(pairs,1)
            i = pairs(k,1);  j = pairs(k,2);
            plot(h, [T(i,1) E(j,1)], [T(i,2), E(j,2)], '-', ...
                'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
        end
    
        plot(h, T(:,1), T(:,2), 'o', 'MarkerSize', 9,  'LineWidth', 1.5, ...
            'Color', [0.85 0.33 0.10],  'DisplayName', 'True');
        plot(h, E(:,1), E(:,2), 'x', 'MarkerSize', 11, 'LineWidth', 1.8, ...
            'Color', [0 0.45 0.74], 'DisplayName', 'Estimated');
    
        if ~isempty(uT)
            plot(h, T(uT,1), T(uT,2), 'o', 'MarkerSize', 15, 'LineWidth', 1.5, ...
                'Color', [0.6 0 0],   'DisplayName', 'Missed');
        end
        if ~isempty(uE)
            plot(h, E(uE,1), E(uE,2), 's', 'MarkerSize', 15, 'LineWidth', 1.5, ...
                'Color', [0.6 0 0.6], 'DisplayName', 'False alarm');
        end
    
        % pair labels at segment midpoints
        for k = 1:size(pairs,1)
            mid = 0.5*(T(pairs(k,1),:) + E(pairs(k,2),:));
            
            txt = sprintf('  %d (%.3g)', k, d(k));
            
            text(h, mid(1), mid(2), txt, 'FontSize', 12, ...
                'VerticalAlignment', 'bottom', 'Interpreter', 'none');
        end
    
        axis(h, 'equal');  grid(h, 'on');  box(h, 'on');
        xlabel('crossrange [m]')
        ylabel('range [m]')
        legend(h, 'Location', 'southwest');
    
        if isempty(d)
            title(h, 'No matched pairs');
        else
            title(h, sprintf('NOMP RMSE %.4g  |  matched %d  |  missed %d  |  FA %d', ...
                sqrt(mean(d.^2)), numel(d), numel(uT), numel(uE)), 'FontSize', font_size);
        end
    
        if ~wasHold
            hold(h, 'off');
        end
        xlim([min(x_array), max(x_array)])
        ylim([min(y_array)+u0, max(y_array)+u0])
        set(gca,'FontSize',font_size)
        axis square
    end

    if isfield(options, 'execute_omp') ...
            && options.execute_omp

        T = target_locations; 
        E = x_hat.omp.positions;
        d = x_hat.omp.d;
        pairs = x_hat.omp.pairs;

        subplot(1,num_plots,iplot)
        iplot = iplot + 1;

        T(:,2) = T(:,2) + u0;
        E(:,2) = E(:,2) + u0;
    
        if ~isempty(pairs)
            if max(pairs(:,1)) > size(T,1) || max(pairs(:,2)) > size(E,1) ...
                    || min(pairs(:)) < 1
                error('pairs contains indices outside the location arrays.');
            end
        end
        
        h = gca;
    
        uT = setdiff((1:size(T,1)).', pairs(:,1));   % unmatched true  -> missed
        uE = setdiff((1:size(E,1)).', pairs(:,2));   % unmatched est.  -> false alarm
    
        wasHold = ishold(h);
        hold(h, 'on');
    
        % connecting segments first so markers draw on top
        for k = 1:size(pairs,1)
            i = pairs(k,1);  j = pairs(k,2);
            plot(h, [T(i,1) E(j,1)], [T(i,2), E(j,2)], '-', ...
                'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
        end
    
        plot(h, T(:,1), T(:,2), 'o', 'MarkerSize', 9,  'LineWidth', 1.5, ...
            'Color', [0.85 0.33 0.10],  'DisplayName', 'True');
        plot(h, E(:,1), E(:,2), 'x', 'MarkerSize', 11, 'LineWidth', 1.8, ...
            'Color', [0 0.45 0.74], 'DisplayName', 'Estimated');
    
        if ~isempty(uT)
            plot(h, T(uT,1), T(uT,2), 'o', 'MarkerSize', 15, 'LineWidth', 1.5, ...
                'Color', [0.6 0 0],   'DisplayName', 'Missed');
        end
        if ~isempty(uE)
            plot(h, E(uE,1), E(uE,2), 's', 'MarkerSize', 15, 'LineWidth', 1.5, ...
                'Color', [0.6 0 0.6], 'DisplayName', 'False alarm');
        end
    
        % pair labels at segment midpoints
        for k = 1:size(pairs,1)
            mid = 0.5*(T(pairs(k,1),:) + E(pairs(k,2),:));
            
            txt = sprintf('  %d (%.3g)', k, d(k));
            
            text(h, mid(1), mid(2), txt, 'FontSize', 12, ...
                'VerticalAlignment', 'bottom', 'Interpreter', 'none');
        end
    
        axis(h, 'equal');  grid(h, 'on');  box(h, 'on');
        xlabel('crossrange [m]')
        ylabel('range [m]')
        legend(h, 'Location', 'southwest');
    
        if isempty(d)
            title(h, 'No matched pairs');
        else
            title(h, sprintf('OMP RMSE %.4g  |  matched %d  |  missed %d  |  FA %d', ...
                sqrt(mean(d.^2)), numel(d), numel(uT), numel(uE)), 'FontSize', font_size);
        end
    
        if ~wasHold
            hold(h, 'off');
        end
        xlim([min(x_array), max(x_array)])
        ylim([min(y_array)+u0, max(y_array)+u0])
        set(gca,'FontSize',font_size)
        axis square
    end

    if isfield(options, 'execute_bp') ...
            && options.execute_bp

        T = target_locations; 
        E = x_hat.bp.positions;
        d = x_hat.bp.d;
        pairs = x_hat.bp.pairs;

        subplot(1,num_plots,iplot)
        iplot = iplot + 1;

        T(:,2) = T(:,2) + u0;
        E(:,2) = E(:,2) + u0;
    
        if ~isempty(pairs)
            if max(pairs(:,1)) > size(T,1) || max(pairs(:,2)) > size(E,1) ...
                    || min(pairs(:)) < 1
                error('pairs contains indices outside the location arrays.');
            end
        end
        
        h = gca;
    
        uT = setdiff((1:size(T,1)).', pairs(:,1));   % unmatched true  -> missed
        uE = setdiff((1:size(E,1)).', pairs(:,2));   % unmatched est.  -> false alarm
    
        wasHold = ishold(h);
        hold(h, 'on');
    
        % connecting segments first so markers draw on top
        for k = 1:size(pairs,1)
            i = pairs(k,1);  j = pairs(k,2);
            plot(h, [T(i,1) E(j,1)], [T(i,2), E(j,2)], '-', ...
                'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
        end
    
        plot(h, T(:,1), T(:,2), 'o', 'MarkerSize', 9,  'LineWidth', 1.5, ...
            'Color', [0.85 0.33 0.10],  'DisplayName', 'True');
        plot(h, E(:,1), E(:,2), 'x', 'MarkerSize', 11, 'LineWidth', 1.8, ...
            'Color', [0 0.45 0.74], 'DisplayName', 'Estimated');
    
        if ~isempty(uT)
            plot(h, T(uT,1), T(uT,2), 'o', 'MarkerSize', 15, 'LineWidth', 1.5, ...
                'Color', [0.6 0 0],   'DisplayName', 'Missed');
        end
        if ~isempty(uE)
            plot(h, E(uE,1), E(uE,2), 's', 'MarkerSize', 15, 'LineWidth', 1.5, ...
                'Color', [0.6 0 0.6], 'DisplayName', 'False alarm');
        end
    
        % pair labels at segment midpoints
        for k = 1:size(pairs,1)
            mid = 0.5*(T(pairs(k,1),:) + E(pairs(k,2),:));
            
            txt = sprintf('  %d (%.3g)', k, d(k));
            
            text(h, mid(1), mid(2), txt, 'FontSize', 12, ...
                'VerticalAlignment', 'bottom', 'Interpreter', 'none');
        end
    
        axis(h, 'equal');  grid(h, 'on');  box(h, 'on');
        xlabel('crossrange [m]')
        ylabel('range [m]')
        legend(h, 'Location', 'southwest');
    
        if isempty(d)
            title(h, 'No matched pairs');
        else
            title(h, sprintf('BP RMSE %.4g  |  matched %d  |  missed %d  |  FA %d', ...
                sqrt(mean(d.^2)), numel(d), numel(uT), numel(uE)), 'FontSize', font_size);
        end
    
        if ~wasHold
            hold(h, 'off');
        end
        xlim([min(x_array), max(x_array)])
        ylim([min(y_array)+u0, max(y_array)+u0])
        set(gca,'FontSize',font_size)
        axis square
    end

    sgtitle(super_title, 'FontSize', 40);
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(f, ['plots/position_', num2str(scenario_number), '.png'])
end

