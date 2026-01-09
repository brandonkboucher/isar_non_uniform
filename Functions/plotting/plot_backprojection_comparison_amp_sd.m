function plot_backprojection_comparison_sd(...
    rx_signal_bp, ...
    x_array, ...
    y_array, ...
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
    
    cmax = max(abs(truth_output.rx_signal_bp.'), [], "all");

    %% original plot
    subplot(1,2,1)
    imagesc(x_array, y_array, abs(truth_output.rx_signal_bp.'))
    colormap("gray")
    xlabel('x (cross-range) [m]', 'FontSize', 16)
    ylabel('y (range) [m]', 'FontSize', 16)
    axis square
    grid on
    c2 = colorbar;
    c2.Label.String = 'Amplitude [Linear]'; 
    c2.FontSize = 16;
    caxis([0, cmax])
    set(gca,'FontSize',16)
    
    % orig_title = create_title('SD image', original_omega, inf);
    % title(orig_title, 'FontSize', 24, 'Interpreter','latex')
    
    %% new plot
    subplot(1,2,2)
    imagesc(x_array, y_array, abs(rx_signal_bp.'));
    colormap("gray")
    xlabel('x (cross-range) [m]', 'FontSize', 16)
    ylabel('y (range) [m]', 'FontSize', 16)
    axis square
    grid on
    c1 = colorbar;
    c1.Label.String = 'Amplitude [Linear]';  
    c1.FontSize = 16;
    caxis([0, cmax])
    set(gca,'FontSize',16)

    % title_str = create_title('SD image', omega, omega_dot);
    % title(title_str, 'FontSize', 24, 'Interpreter','latex')

    set(gcf, 'Position', get(0, 'Screensize'));
    if save_plots
        saveas(f, ['plots/backproj_sd_amp_comp_', num2str(omega_dot), '.png'])
    end
end

