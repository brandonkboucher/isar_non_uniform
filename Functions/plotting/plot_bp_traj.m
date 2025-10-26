function plot_bp_traj( ...
    grid_extent, ... % [cross range extent, range extent]
    yaws, ...
    output_struct)

    % extract target positions, scatterer positions and the
    % output image
    target_positions = output_struct.target_positions;
    scatterer_positions = output_struct.scatterer_positions;
    image = output_struct.rx_bp;
    ranges = output_struct.ranges;
    range_target2 = ranges(:,2);

    buffer = norm(grid_extent/2);
    ymin = target_positions(1,2) - buffer;
    ymax = target_positions(1,2) + buffer;
    xmin = target_positions(1,1) - buffer;
    xmax = target_positions(1,1) + buffer;

    grid_vertices = ...
        [grid_extent(1)/2, grid_extent(2)/2, 0; ... % quad 1
         grid_extent(1)/2, -grid_extent(2)/2, 0; ... % quad 2
         -grid_extent(1)/2, -grid_extent(2)/2, 0; ... % quad 3
         -grid_extent(1)/2, grid_extent(2)/2, 0; ... % quad 4
         grid_extent(1)/2, grid_extent(2)/2, 0]; % close square


    j = 1;
    f = figure('Visible','off');

    subplot(1,2,1)

    imagesc(...
        output_struct.x_bp + target_positions(1,1), ...
        output_struct.y_bp + target_positions(1,2), ...
        20*log10(abs(output_struct.rx_bp.') + eps));

    title('Backprojection image', 'FontSize', 24)
    xlabel('x (cross-range)', 'FontSize', 16)
    ylabel('y (range)', 'FontSize', 16)
    axis square
    colorbar   
    set(gca,'FontSize',16)

    subplot(1,2,2)
    for ipulse = 1:size(yaws, 1)-1

        % plot the Backprojection grid vertices to ensure
        % they line up with the rotation of the target
        yaw = yaws(ipulse);
        rot = ...
            [cos(yaw), -sin(yaw), 0; ...
             sin(yaw), cos(yaw), 0; ...
             0, 0, 1];

        rotated_grid_vertices = (rot * grid_vertices')';
        rotated_grid_vertices = rotated_grid_vertices ...
            + target_positions(ipulse, :);
        
        g = plot(rotated_grid_vertices(:,1)', rotated_grid_vertices(:,2)', '-k');
        hold on
        
        scatterers = squeeze(scatterer_positions(ipulse, :, :));

        s = plot(scatterers(:,1)', scatterers(:,2)', 'Marker','o', 'LineStyle','none', 'LineWidth', 5);
        
        r = plot([0,scatterers(2,1)], [0,scatterers(2,2)], 'LineWidth',2);
        
        t = text(scatterers(2,1)+1, scatterers(2,2), 0, sprintf('range = %4.2f', range_target2(ipulse)));
        
        ylim([ymin, ymax])
        xlim([xmin, xmax])
        
        set(gcf, 'Position', get(0, 'Screensize'));
        set(gca, 'YDir', 'reverse')
        F(j) = getframe(gcf);
        j = j + 1;
        delete(g)
        delete(s)
        delete(r)
        delete(t)
    end

    % create the video writer with 1 fps
    writerObj = VideoWriter('bp_target_trajectory', 'MPEG-4');

    % open the video writer
    open(writerObj);
    % write the frames to the video
    for i=1:length(F)
        % convert the image to a frame
        frame = F(i) ;    
        writeVideo(writerObj, frame);
    end
    % close the writer object
    close(writerObj);

    % xmin = min([scatterer_positions(:,:,1)], [], "all");
    % ymin = min([scatterer_positions(:,:,2)], [], "all");
    % 
    % if plotting_settings.visible
    %     figure
    % else
    %     f = figure('Visible','off');
    % end
    % 
    % for ipt = 1:size(scatterer_positions, 2)
    % 
    %     plot( ...
    %         scatterer_positions(:,ipt,1), ...
    %         scatterer_positions(:,ipt,2), ...
    %         'DisplayName', 'Target Trajectory', 'Marker','+', 'LineStyle','none')
    % 
    %     hold on
    % 
    % 
    % end
    % 
    % xlim([0, max(scatterer_positions(:,:,1),[], "all")])
    % ylim([0, max(scatterer_positions(:,:,2),[], "all")])
    % 
    % title('Target Trajectory', 'FontSize', plotting_settings.title_font_size)
    % xlabel('X [m]', 'FontSize', plotting_settings.label_font_size)
    % ylabel('Y [m]', 'FontSize', plotting_settings.label_font_size)
    % grid on
    % axis square
    % set(gca,'FontSize',plotting_settings.axis_font_size)
    % 
    % if ~obj.visible
    %     set(gcf, 'Position', get(0, 'Screensize'));
    %     saveas(f, [plotting_settings.plot_folder, 'trajectory2d.png'])
    % end

end

