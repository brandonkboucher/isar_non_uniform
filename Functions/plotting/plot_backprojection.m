function plot_backprojection(...
    rx_signal_bp, ...
    x_array, ...
    y_array, ...
    show_plots, ...
    save_plots)
    
    visible = 'off';
    if show_plots
        visible = 'on'; 
    end

    f = figure('Visible',visible);
    subplot(1,2,1)
    imagesc(x_array, y_array, 20*log10(abs(rx_signal_bp.') + eps));
    xlabel('x (cross-range)', 'FontSize', 16)
    ylabel('y (range)', 'FontSize', 16)
    axis square
    grid on
    colorbar   
    set(gca,'FontSize',16)
    title('Backprojection image - Log scaled', 'FontSize', 24)
    
    
    subplot(1,2,2)
    imagesc(x_array, y_array, abs(rx_signal_bp.'))
    xlabel('x (cross-range)', 'FontSize', 16)
    ylabel('y (range)', 'FontSize', 16)
    axis square
    grid on
    colorbar   
    set(gca,'FontSize',16)
    title('Backprojection image', 'FontSize', 24)
    
    set(gcf, 'Position', get(0, 'Screensize'));
    if save_plots
        saveas(f, 'plots/backproj.png')
    end
end

