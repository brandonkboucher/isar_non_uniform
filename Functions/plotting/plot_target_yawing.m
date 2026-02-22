function plot_target_yawing(...
    t_slow, ...
    target, ...
    show_plots, ...
    save_plots, ...
    include_title, ...
    plotting_dir)
    
    % range and trajectory
    visible = 'off';
    if show_plots
        visible = 'on'; 
    end
    
    N = 3;
    if isprop(target, 'yawing_jerk')
        N = 4;
    end

    f = figure('Visible',visible);
    subplot(1,N,1)
    plot(t_slow, target.yaws, 'LineWidth', 2)
    
    xlabel('Slow time [s]', 'FontSize', 16)
    ylabel('Yaw [rad]', 'FontSize', 16)
    grid on
    axis square
    ax = gca;
    set(ax,'FontSize',16)
    if include_title
        title('Target Yaws', 'FontSize', 24)
    end

    subplot(1,N,2)
    plot(t_slow, target.yawing_rate, 'LineWidth', 2)
    
    xlabel('Slow time [s]', 'FontSize', 16)
    ylabel('Yawing Rate [rad/s]', 'FontSize', 16)
    grid on
    axis square
    ax = gca;
    set(ax,'FontSize',16)
    if include_title
        title('Target Yawing Rate', 'FontSize', 24)
    end

    subplot(1,N,3)
    plot(t_slow, target.yawing_acceleration, 'LineWidth', 2)
    
    xlabel('Slow time [s]', 'FontSize', 16)
    ylabel('Yawing Acceleration [rad/s/s]', 'FontSize', 16)
    grid on
    axis square
    ax = gca;
    set(ax,'FontSize',16)
    if include_title
        title('Target Yawing Acceleration', 'FontSize', 24)
    end

    if N > 3
        subplot(1,N,4)
        plot(t_slow, target.yawing_jerk, 'LineWidth', 2)
        
        xlabel('Slow time [s]', 'FontSize', 16)
        ylabel('Yawing Jerk [rad/s/s/s]', 'FontSize', 16)
        grid on
        axis square
        ax = gca;
        set(ax,'FontSize',16)
        if include_title
            title('Target Yawing Jerk', 'FontSize', 24)
        end
    end

    set(gcf, 'Position', get(0, 'Screensize'));
    if save_plots
        saveas(f,[plotting_dir,'/target_rotation.png'])
    end
end

