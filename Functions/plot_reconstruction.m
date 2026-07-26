
function plot_reconstruction(...
    scenario_number, ...
    target_locations,...
    u0,...
    x_array, ...
    y_array, ...
    x_hat, ...
    options, ...
    super_title)

    num_plots = ...
        options.execute_itsa ...
        + options.execute_omp ...
        + options.execute_sbl ...
        + options.execute_fbp;

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

    if options.execute_omp
        subplot(1,num_plots,iplot)
        iplot = iplot + 1;
        imagesc(x_array, y_array + u0, abs(x_hat.omp))
        axis square
        hold on; overlay_truth(); hold off
        xlabel('crossrange [m]')
        ylabel('range [m]')
        title('OMP reconstruction (o = true off-grid positions)', 'FontSize', 24)
        set(gca, "FontSize",16)
        colormap gray
        colorbar
    end
    
    if options.execute_itsa
        subplot(1,num_plots,iplot)
        iplot = iplot + 1;
        imagesc(x_array, y_array + u0, abs(x_hat.lasso))
        axis square
        hold on; overlay_truth(); hold off
        xlabel('crossrange [m]')
        ylabel('range [m]')
        title('LASSO (ISTA) reconstruction (o = true off-grid positions)', 'FontSize', 24)
        set(gca, "FontSize",16)
        colormap gray
        colorbar
    end
    
    if options.execute_sbl
        subplot(1,num_plots,iplot)
        iplot = iplot + 1;
        imagesc(x_array, y_array + u0, abs(x_hat.sbl))
        axis square
        hold on; overlay_truth(); hold off
        xlabel('crossrange [m]')
        ylabel('range [m]')
        title('EM-SBL reconstruction (o = true off-grid positions)', 'FontSize', 24)
        set(gca, "FontSize",16)
        colormap gray
        colorbar
    end
    
    if options.execute_fbp
        subplot(1,num_plots,iplot)
        iplot = iplot + 1;
        imagesc(x_array, y_array + u0, abs(x_hat.fbp))
        axis square
        hold on; overlay_truth(); hold off
        xlabel('crossrange [m]')
        ylabel('range [m]')
        title('BP reconstruction (o = true off-grid positions)', 'FontSize', 24)
        set(gca, "FontSize",16)
        colormap gray
        colorbar
    end
    sgtitle(super_title, 'FontSize', 40);
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(f, ['plots/scenario_', num2str(scenario_number), '.png'])
end

