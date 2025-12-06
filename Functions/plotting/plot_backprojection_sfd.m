function plot_backprojection_sfd(...
    rx_signal_sfd, ...
    kx, ...
    ky, ...
    show_plots, ...
    save_plots, ...
    filename)
    
    visible = 'off';
    if show_plots
        visible = 'on'; 
    end

    f = figure('Visible',visible);
    subplot(1,2,1)
    imagesc(kx, ky, 20*log10(abs(rx_signal_sfd.') + eps));
    colormap("gray")
    colorbar
    xlabel('kx [1/m]', 'FontSize', 16)
    ylabel('ky [1/m]', 'FontSize', 16)
    axis square
    grid on
    c1 = colorbar;
    c1.Label.String = 'Amplitude [dB]';   
    set(gca,'FontSize',16)
    title('Backprojection SFD image - Log scaled', 'FontSize', 24)
    
    
    subplot(1,2,2)
    imagesc(kx, ky, abs(rx_signal_sfd.'))
    colormap("gray")
    colorbar
    xlabel('kx [1/m]', 'FontSize', 16)
    ylabel('ky [1/m]', 'FontSize', 16)
    axis square
    grid on
    c2 = colorbar;
    c2.Label.String = 'Amplitude [Linear]';   
    set(gca,'FontSize',16)
    title('Backprojection SFD image', 'FontSize', 24)
    
    set(gcf, 'Position', get(0, 'Screensize'));
    if save_plots
        saveas(f, ['plots/backproj_sfd_', filename, '.png'])
    end
end

