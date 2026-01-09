function plot_target_range_and_traj(...
    t_slow, ...
    ranges, ...
    scatterer_positions, ...
    yaws, ...
    show_plots, ...
    save_plots, ...
    filename)
    
    % range and trajectory
    visible = 'off';
    if show_plots
        visible = 'on'; 
    end
    
    f = figure('Visible',visible);
    subplot(1,3,1)
    plot(t_slow, ranges, 'LineWidth', 2)
    xlabel('Slow time', 'FontSize', 16)
    ylabel('Range', 'FontSize', 16)
    axis square
    grid on
    ax = gca;
    set(ax,'FontSize',16)
    %title('Target Range', 'FontSize', 24)
    ax.YDir = "reverse";
    
    subplot(1,3,2)
    for ipt = 1:size(scatterer_positions, 2)
    
        plot( ...
            scatterer_positions(:,ipt,1), ...
            scatterer_positions(:,ipt,2), ...
            'DisplayName', 'Target Trajectory', ...
            'Marker','+', 'LineStyle','none')
        
        hold on
    
    end
    
    ylim([-10, 10])
    xlim([-10, 10])
    
    xlabel('X [m]', 'FontSize', 16)
    ylabel('Y [m]', 'FontSize', 16)
    grid on
    axis square
    ax = gca;
    set(ax,'FontSize',16)
    %title('Target Trajectory', 'FontSize', 24)
    
    ax.YDir = "reverse";

    subplot(1,3,3)
    plot(t_slow, yaws, 'LineWidth', 2)
    
    xlabel('Slow time [s]', 'FontSize', 16)
    ylabel('Target Yaw [rad]', 'FontSize', 16)
    grid on
    axis square
    ax = gca;
    set(ax,'FontSize',16)
    %title('Target Rotation', 'FontSize', 24)

    set(gcf, 'Position', get(0, 'Screensize'));
    if save_plots
        saveas(f,['plots/range_and_trj', filename, '.png'])
    end
end

