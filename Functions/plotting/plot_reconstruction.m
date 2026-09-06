
function plot_reconstruction(...
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

    if isfield(options, 'execute_mod_omp') ...
            && options.execute_mod_omp
        num_plots = num_plots + 1;
    end

    if isfield(options, 'execute_nomp') ...
            && options.execute_nomp
        num_plots = num_plots + 1;
    end

    if isfield(options, 'execute_promp') ...
            && options.execute_promp
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

    % panel titles crowd each other as algorithms are added to the row, so
    % scale them down once the row holds more than three panels
    title_font_size = font_size;
    if num_plots > 3
        title_font_size = round(font_size * 3 / num_plots);
    end

    % true (off-grid) scatterer positions in plot coordinates: crossrange on x,
    % range (offset by u0) on y. These continuous positions do not lie on grid
    % points, so they are overlaid as markers rather than shown as an image.
    true_x = target_locations(:,1);
    true_y = target_locations(:,2) + u0;

    % helper to overlay the true positions on the current axes
    overlay_truth = @() plot(true_x, true_y, 'o', ...
        'MarkerEdgeColor', [1 0 0], 'MarkerSize', 12, 'LineWidth', 1.5);
    
    f = figure('Visible','off');
    iplot = 1;

    if isfield(options, 'execute_nomp') ...
            && options.execute_nomp

        subplot(1,num_plots,iplot)
        iplot = iplot + 1;

        x_nomp = x_hat.nomp.positions(:,1);
        y_nomp = x_hat.nomp.positions(:,2);
        alpha = x_hat.nomp.alpha;

        % helper to overlay the true positions on the current axes
        plot(x_nomp, y_nomp + u0, '*', ...
            'MarkerEdgeColor', [0 0 1], ...
            'MarkerSize', 12, ...
            'LineWidth', 1.5);

        axis square
        hold on; overlay_truth(); hold off
        xlabel('crossrange [m]')
        ylabel('range [m]')
        title(sprintf('NOMP reconstruction (error: %.2e)', x_hat.nomp.error), 'FontSize', title_font_size)
        set(gca, "FontSize",font_size)
        xlim([min(x_array), max(x_array)])
        ylim([min(y_array)+u0, max(y_array)+u0])
        set(gca, 'YDir', 'reverse');
        set(gca,'FontSize',font_size)
    end

    if isfield(options, 'execute_promp') ...
            && options.execute_promp

        subplot(1,num_plots,iplot)
        iplot = iplot + 1;

        x_promp = x_hat.promp.positions(:,1);
        y_promp = x_hat.promp.positions(:,2);
        alpha = x_hat.promp.alpha;

        % helper to overlay the true positions on the current axes
        plot(x_promp, y_promp + u0, '*', ...
            'MarkerEdgeColor', [0 0.5 0], ...
            'MarkerSize', 12, ...
            'LineWidth', 1.5);

        axis square
        hold on; overlay_truth(); hold off
        xlabel('crossrange [m]')
        ylabel('range [m]')
        title(sprintf('PROMP reconstruction (error: %.2e)', x_hat.promp.error), 'FontSize', title_font_size)
        set(gca, "FontSize",font_size)
        xlim([min(x_array), max(x_array)])
        ylim([min(y_array)+u0, max(y_array)+u0])
        set(gca, 'YDir', 'reverse');
        set(gca,'FontSize',font_size)
    end

    if isfield(options, 'execute_omp') ...
            && options.execute_omp

        subplot(1,num_plots,iplot)
        iplot = iplot + 1;
        if options.log_scale_plotting
            imagesc(x_array, y_array + u0, 20*log10(abs(x_hat.omp.image)))
        else
            imagesc(x_array, y_array + u0, abs(x_hat.omp.image))
        end
        axis square
        hold on; overlay_truth(); hold off
        xlabel('crossrange [m]')
        ylabel('range [m]')
        title(sprintf('OMP reconstruction (error: %.2e)', x_hat.omp.error), 'FontSize', title_font_size)
        set(gca, "FontSize",font_size)
        colormap gray
        c2 = colorbar;
        if options.log_scale_plotting
            c2.Label.String = 'Log-Scaled [dB]';
        else
            c2.Label.String = 'Amplitude [Linear]'; 
        end
        c2.FontSize = font_size;
        set(gca,'FontSize',font_size)
    end
    
    if isfield(options, 'execute_mod_omp') ...
            && options.execute_mod_omp

        subplot(1,num_plots,iplot)
        iplot = iplot + 1;

        % the modified OMP image is [Ny, num_of_amb * Nx]: one grid image per
        % ambiguity stacked along crossrange. It therefore carries its own
        % crossrange axis, which spans several unambiguous extents Wx rather
        % than the single extent of `x_array`.
        x_array_mod = x_hat.mod_omp.x_array;
        amb_index   = x_hat.mod_omp.amb_index;
        Nx_block    = size(x_hat.mod_omp.image,2) / numel(amb_index);

        if options.log_scale_plotting
            imagesc(x_array_mod, y_array + u0, 20*log10(abs(x_hat.mod_omp.image)))
        else
            imagesc(x_array_mod, y_array + u0, abs(x_hat.mod_omp.image))
        end

        % keep the same panel box as the single-ambiguity images so the row
        % stays uniform; the stacked image is simply compressed in crossrange
        axis square
        hold on
        overlay_truth();

        % separate the ambiguity blocks and label each with its index. The
        % y-axis is reversed by imagesc, so min(y_array) is the top edge.
        for iamb = 1:numel(amb_index)
            if iamb > 1
                xb = mean(x_array_mod(round((iamb-1)*Nx_block) + [0 1]));
                plot([xb xb], [min(y_array) max(y_array)] + u0, '--', ...
                    'Color', [1 1 0], 'LineWidth', 1.5);
            end
            xc = mean(x_array_mod(round((iamb-1)*Nx_block) + [1 Nx_block]));
            text(xc, min(y_array) + u0, sprintf('%+d', amb_index(iamb)), ...
                'Color', [1 1 0], 'FontSize', font_size, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
        end
        hold off

        xlabel('crossrange [m]')
        ylabel('range [m]')
        title(sprintf('Modified OMP reconstruction (error: %.2e)', ...
            x_hat.mod_omp.error), 'FontSize', title_font_size)
        set(gca, "FontSize",font_size)
        xlim([min(x_array_mod), max(x_array_mod)])
        ylim([min(y_array)+u0, max(y_array)+u0])
        colormap gray
        c2 = colorbar;
        if options.log_scale_plotting
            c2.Label.String = 'Log-Scaled [dB]';
        else
            c2.Label.String = 'Amplitude [Linear]'; 
        end
        c2.FontSize = font_size;
        set(gca,'FontSize',font_size)
    end
    
    if isfield(options, 'execute_itsa') ...
            && options.execute_itsa
        subplot(1,num_plots,iplot)
        iplot = iplot + 1;
        if options.log_scale_plotting
            imagesc(x_array, y_array + u0, 20*log10(abs(x_hat.lasso.image)))
        else
            imagesc(x_array, y_array + u0, abs(x_hat.lasso.image))
        end
        axis square
        hold on; overlay_truth(); hold off
        xlabel('crossrange [m]')
        ylabel('range [m]')
        title(sprintf('LASSO (ISTA) reconstruction (error: %.2e)', x_hat.lasso.error), 'FontSize', title_font_size)
        set(gca, "FontSize",font_size)
        colormap gray
        c2 = colorbar;
        if options.log_scale_plotting
            c2.Label.String = 'Log-Scaled [dB]';
        else
            c2.Label.String = 'Amplitude [Linear]'; 
        end
        c2.FontSize = font_size;
        set(gca,'FontSize',font_size)
    end
    
    if isfield(options, 'execute_sbl') ...
            && options.execute_sbl
        subplot(1,num_plots,iplot)
        iplot = iplot + 1;
        if options.log_scale_plotting
            imagesc(x_array, y_array + u0, 20*log10(abs(x_hat.sbl.image)))
        else
            imagesc(x_array, y_array + u0, abs(x_hat.sbl.image))
        end
        axis square
        hold on; overlay_truth(); hold off
        xlabel('crossrange [m]')
        ylabel('range [m]')
        title(sprintf('EM-SBL reconstruction (error: %.2e)', x_hat.sbl.error), 'FontSize', title_font_size)
        set(gca, "FontSize",font_size)
        colormap gray
        c2 = colorbar;
        if options.log_scale_plotting
            c2.Label.String = 'Log-Scaled [dB]';
        else
            c2.Label.String = 'Amplitude [Linear]'; 
        end 
        c2.FontSize = font_size;
        set(gca,'FontSize',font_size)
    end
    
    if isfield(options, 'execute_bp') ...
            && options.execute_bp

        subplot(1,num_plots,iplot)
        iplot = iplot + 1;
        if options.log_scale_plotting
            imagesc(x_array, y_array + u0, 20*log10(abs(x_hat.bp.image)))
        else
            imagesc(x_array, y_array + u0, abs(x_hat.bp.image))
        end
        axis square
        hold on; overlay_truth(); hold off
        xlabel('crossrange [m]')
        ylabel('range [m]')
        title(sprintf('BP normalized reconstruction (error: %.2e)', x_hat.bp.error), 'FontSize', title_font_size)
        set(gca, "FontSize",font_size)
        colormap gray
        c2 = colorbar;
        if options.log_scale_plotting
            c2.Label.String = 'Log-Scaled [dB]';
        else
            c2.Label.String = 'Amplitude [Linear]'; 
        end
        c2.FontSize = font_size;
        set(gca,'FontSize',font_size)
    end
    % the per-panel `set(gca,'FontSize',...)` calls above reset the titles, so
    % apply the scaled title size once the whole row has been drawn
    ax = findall(f, 'Type', 'axes');
    for iax = 1:numel(ax)
        if ~isempty(ax(iax).Title)
            set(ax(iax).Title, 'FontSize', title_font_size);
        end
    end

    sgtitle(super_title, 'FontSize', 40);
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(f, ['plots/scenario_', num2str(scenario_number), '.png'])
end

