function plot_backprojection_sd(...
    rx_signal_bp, ...
    x_array, ...
    y_array, ...
    show_plots, ...
    save_plots, ...
    varargin)
    
    visible = 'off';
    if show_plots
        visible = 'on'; 
    end

    f = figure('Visible',visible);
    subplot(1,2,1)
    imagesc(x_array, y_array, 20*log10(abs(rx_signal_bp.') + eps));
    colormap("gray")
    
    xlabel('x (cross-range)', 'FontSize', 16)
    ylabel('y (range)', 'FontSize', 16)
    axis square
    grid on
    c1 = colorbar;
    c1.Label.String = 'Amplitude [dB]';
    set(gca,'FontSize',16)
    %title('SD image - Log scaled', 'FontSize', 24)
    
    
    subplot(1,2,2)
    imagesc(x_array, y_array, abs(rx_signal_bp.'))
    colormap("gray")
    colorbar
    xlabel('x (cross-range)', 'FontSize', 16)
    ylabel('y (range)', 'FontSize', 16)
    axis square
    grid on
    c2 = colorbar;
    c2.Label.String = 'Amplitude [Linear]';    
    set(gca,'FontSize',16)
    %title('SD image', 'FontSize', 24)
    
    set(gcf, 'Position', get(0, 'Screensize'));
    if save_plots
        
        saveas(f, 'plots/backproj_sd.png')

    end
end

