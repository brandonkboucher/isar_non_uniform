function plot_backproj_comparison_amp_sfd(...
    rx_signal_sfd, ...
    kx, ...
    ky, ...
    show_plots, ...
    save_plots, ...
    filename, ...
    truth_data_dir)
    
    visible = 'off';
    if show_plots
        visible = 'on'; 
    end

    % load the truth data FOR omega_dot = pi/8
    truth_output = load(truth_data_dir);
    truth_output = truth_output.output;
    
    cmax = max(abs(truth_output.rx_signal_sfd.'), [], "all");

    f = figure('Visible',visible);
    subplot(1,2,1)
    imagesc(kx, ky, abs(truth_output.rx_signal_sfd.'))
    colormap("gray")
    xlabel('kx [1/m]', 'FontSize', 16)
    ylabel('ky [1/m]', 'FontSize', 16)
    axis square
    grid on
    c2 = colorbar;
    c2.Label.String = 'Amplitude [Linear]'; 
    c2.FontSize = 16;
    caxis([0, cmax])
    set(gca,'FontSize',16)
    
    title('Backprojection SFD image ($\omega = \pi / 4$ rad/s, $\dot{\omega} = 0$ rad/s/s)', 'FontSize', 24, 'Interpreter','latex')

    
    subplot(1,2,2)
    imagesc(kx, ky, abs(rx_signal_sfd.'))
    colormap("gray")
    xlabel('kx [1/m]', 'FontSize', 16)
    ylabel('ky [1/m]', 'FontSize', 16)
    axis square
    grid on
    c1 = colorbar;
    c1.Label.String = 'Amplitude [Linear]';
    c1.FontSize = 16;
    caxis([0, cmax])
    set(gca,'FontSize',16)

    factor = str2double(filename);
    if factor > 1 && factor ~= inf

        title_str = ['Backprojection SFD image ($\omega = \pi / 4$ rad/s, $\dot{\omega} = \pi/', num2str(factor),'$ rad/s/s)'];
    elseif factor == inf
        title_str = 'Backprojection SFD image ($\omega = \pi / 4$ rad/s, $\dot{\omega} = 0$ rad/s/s)';
    else
        factor = 1 / factor;
        title_str = ['Backprojection SFD image ($\omega = \pi / 4$ rad/s, $\dot{\omega} = ', num2str(factor), '\pi$ rad/s/s)'];
    end
    title(title_str, 'FontSize', 24, 'Interpreter','latex')


    set(gcf, 'Position', get(0, 'Screensize'));
    if save_plots
        saveas(f, ['plots/backproj_sfd_comp_amp_', filename, '.png'])
    end
end

