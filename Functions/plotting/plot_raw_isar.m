function plot_raw_isar(...
    rx_signal, ...
    range_array, ...
    show_plots, ...
    save_plots)
    
    visible = 'off';
    if show_plots
        visible = 'on'; 
    end

    f = figure('Visible',visible);
    subplot(1,2,1)
    imagesc(1:size(rx_signal,1),range_array, 20*log10(abs(rx_signal.') + eps));
    colormap("gray")
    colorbar
    xlabel('pulse index', 'FontSize', 16)
    ylabel('range [m]', 'FontSize', 16)
    axis square
    grid on
    colorbar   
    set(gca,'FontSize',16)
    title('Raw ISAR image - Log scaled', 'FontSize', 24)
    
    
    subplot(1,2,2)
    imagesc(1:size(range_array,1), range_array, abs(rx_signal.'))
    colormap("gray")
    colorbar
    xlabel('pulse index', 'FontSize', 16)
    ylabel('range [m]', 'FontSize', 16)
    axis square
    grid on
    colorbar   
    set(gca,'FontSize',16)
    title('Raw ISAR image', 'FontSize', 24)
    
    set(gcf, 'Position', get(0, 'Screensize'));
    if save_plots
        saveas(f, 'plots/raw_isar.png')
    end
end

