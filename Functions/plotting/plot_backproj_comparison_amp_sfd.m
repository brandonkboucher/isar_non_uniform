function plot_backproj_comparison_amp_sfd(...
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

    % load the truth data FOR omega_dot = pi/8
    truth_output = load(truth_data_dir);
    truth_output = truth_output.output;
    
    cmax = max(abs(truth_output.rx_signal_sfd.'), [], "all");

    %% original
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
    
    orig_title = create_title('SFD amplitude', original_omega, inf);
    title(orig_title, 'FontSize', 24, 'Interpreter','latex')
    
    %% new
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

    title_str = create_title('SFD amplitude', omega, omega_dot);
    title(title_str, 'FontSize', 24, 'Interpreter','latex')

    set(gcf, 'Position', get(0, 'Screensize'));
    if save_plots
        saveas(f, ['plots/backproj_sfd_comp_amp_', num2str(omega_dot), '.png'])
    end
end

