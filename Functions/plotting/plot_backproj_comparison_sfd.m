function plot_backproj_comparison_sfd(...
    rx_signal_sfd, ...
    kx, ...
    ky, ...
    show_plots, ...
    save_plots, ...
    omega, ...
    omega_dot, ...
    original_omega, ...
    truth_data_dir)
    
    visible = 'off';
    if show_plots
        visible = 'on'; 
    end

    f = figure('Visible',visible);

    % load the truth data FOR omega_dot = pi/8
    truth_output = load(truth_data_dir);
    truth_output = truth_output.output;
    
    cmax = max(20*log10(abs(truth_output.rx_signal_sfd.')), [], "all");

    %% original
    subplot(1,2,1)
    imagesc(kx, ky, 20*log10(abs(truth_output.rx_signal_sfd.')))
    colormap("gray")
    xlabel('kx [1/m]', 'FontSize', 16)
    ylabel('ky [1/m]', 'FontSize', 16)
    axis square
    grid on
    c2 = colorbar;
    c2.Label.String = 'Amplitude [dB]'; 
    c2.FontSize = 16;
    caxis([0, cmax])
    set(gca,'FontSize',16)

    orig_title = create_title('SFD image', original_omega, inf);
    title(orig_title, 'FontSize', 24, 'Interpreter','latex')
    
    %% new
    subplot(1,2,2)
    imagesc(kx, ky, 20*log10(abs(rx_signal_sfd.')))
    colormap("gray")
    xlabel('kx [1/m]', 'FontSize', 16)
    ylabel('ky [1/m]', 'FontSize', 16)
    axis square
    grid on
    c1 = colorbar;
    c1.Label.String = 'Amplitude [dB]'; 
    c1.FontSize = 16;
    caxis([0, cmax])
    set(gca,'FontSize',16)

    title_str = create_title('SFD image', omega, omega_dot);
    title(title_str, 'FontSize', 24, 'Interpreter','latex')

    set(gcf, 'Position', get(0, 'Screensize'));
    if save_plots
        saveas(f, ['plots/backproj_sfd_comp_', num2str(omega_dot), '.png'])
    end
end

