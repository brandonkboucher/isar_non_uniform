function plot_psf(...
    rx_signal_bp, ...
    x_array, ...
    y_array, ...
    show_plots, ...
    save_plots, ...
    omega, ...
    omega_dot)

    rx_signal_bp_norm = ...
        rx_signal_bp / max((rx_signal_bp), [], "all");

    visible = 'off';
    if show_plots
        visible = 'on'; 
    end

    f = figure('Visible',visible);
    subplot(1,2,1)
    surf(x_array, y_array, abs(rx_signal_bp_norm))
    xlabel('x (cross-range)', 'FontSize', 16)
    ylabel('y (range)', 'FontSize', 16)
    zlabel('Normalized Magnitude', 'FontSize', 16)
    axis square
    grid on
    %c1 = colorbar;
    %c1.Label.String = 'Normalized Magnitude';
    %ax = gca;
    %ax.Position = ax.Position - [0 0 .1 .1];
    %c1.Position = c1.Position - [.125 0 0 0];

    set(gca,'FontSize',16)
    
    title_str = create_title('PSF', omega, omega_dot);


    title(title_str, 'FontSize', 24, 'Interpreter','latex')

    subplot(1,2,2)
    imagesc(x_array, y_array, abs(rx_signal_bp_norm))
    xlabel('x (cross-range)', 'FontSize', 16)
    ylabel('y (range)', 'FontSize', 16)
    axis square
    grid on
    c1 = colorbar;
    c1.Label.String = 'Normalized Magnitude';
    set(gca,'FontSize',16)
    title(title_str, 'FontSize', 24, 'Interpreter','latex')
    set(gcf, 'Position', get(0, 'Screensize'));

    if save_plots
        saveas(f, ['plots/psf_', num2str(omega_dot), '.png'])
    end
    
end

