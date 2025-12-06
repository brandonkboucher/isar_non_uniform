function plot_backprojection_comparison_sd(...
    rx_signal_bp, ...
    x_array, ...
    y_array, ...
    show_plots, ...
    save_plots, ...
    filename, ...
    truth_data_dir)
    
    visible = 'off';
    if show_plots
        visible = 'on'; 
    end

    f = figure('Visible',visible);

    % load the truth data FOR omega_dot = pi/8
    truth_output = load(truth_data_dir);
    truth_output = truth_output.output;
    
    cmax = max(20*log10(abs(truth_output.rx_signal_bp.') + eps), [], "all");

    subplot(1,2,1)
    imagesc(x_array, y_array, 20*log10(abs(truth_output.rx_signal_bp.') + eps))
    colormap("gray")
    xlabel('x (cross-range) [m]', 'FontSize', 16)
    ylabel('y (range) [m]', 'FontSize', 16)
    axis square
    grid on
    c2 = colorbar;
    c2.Label.String = 'Amplitude [dB]'; 
    c2.FontSize = 16;
    caxis([0, cmax])
    set(gca,'FontSize',16)
    
    title('Backprojection image ($\omega = \pi / 4$ rad/s, $\dot{\omega} = 0$ rad/s/s)', 'FontSize', 24, 'Interpreter','latex')
    


    subplot(1,2,2)
    imagesc(x_array, y_array, 20*log10(abs(rx_signal_bp.') + eps));
    colormap("gray")
    xlabel('x (cross-range) [m]', 'FontSize', 16)
    ylabel('y (range) [m]', 'FontSize', 16)
    axis square
    grid on
    c1 = colorbar;
    c1.Label.String = 'Amplitude [dB]';   
    c1.FontSize = 16;
    caxis([0, cmax])
    set(gca,'FontSize',16)

    factor = str2double(filename);
    if factor > 1 && factor ~= inf

        title_str = ['Backprojection image ($\omega = \pi / 4$ rad/s, $\dot{\omega} = \pi/', num2str(factor),'$ rad/s/s)'];
    elseif factor == inf
        title_str = 'Backprojection image ($\omega = \pi / 4$ rad/s, $\dot{\omega} = 0$ rad/s/s)';
    else
        factor = 1 / factor;
        title_str = ['Backprojection image ($\omega = \pi / 4$ rad/s, $\dot{\omega} = ', num2str(factor), '\pi$ rad/s/s)'];
    end
    title(title_str, 'FontSize', 24, 'Interpreter','latex')

   
    set(gcf, 'Position', get(0, 'Screensize'));
    if save_plots
        saveas(f, ['plots/backproj_sd_comp_', filename, '.png'])
    end
end

