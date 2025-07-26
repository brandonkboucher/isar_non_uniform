classdef plotting < handle
    %PLOTTING Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        visible = false
        bool_plot_all = false

        bool_plot_range = false
        bool_plot_target_trajectory = false
        bool_video_target_trajectory = false
        bool_plot_target = false

        bool_plot_raw_isar = false
        bool_plot_range_compressed = false
        bool_plot_range_tracking = false
        bool_plot_autofocus = false
        bool_plot_rd = false

        range_array

        title_font_size = 24;
        label_font_size = 16;

    end
    
    methods

        function obj = plotting()
            
        end

        function plot(obj, output_struct)

            if obj.bool_plot_range ...
                || obj.bool_plot_target_trajectory ...
                || obj.bool_video_target_trajectory ...
                || obj.bool_plot_raw_isar ...
                || obj.bool_plot_range_compressed ...
                || obj.bool_plot_range_tracking ...
                || obj.bool_plot_autofocus ...
                || obj.bool_plot_rd ...
                || obj.bool_plot_target ...
                || obj.bool_plot_all

                delete plots/*.png
            end

            if isfield(output_struct, 'range_array')
                obj.range_array = output_struct.range_array;
            end

            if (obj.bool_plot_target_trajectory || obj.bool_plot_all) ...
                    && isfield(output_struct, 'scatterer_positions')
                obj.plot_trajectory(output_struct.scatterer_positions)
            end

            if (obj.bool_video_target_trajectory)...
                && isfield(output_struct, 'scatterer_positions') ...
                && isfield(output_struct, 't_m') ...
                && isfield(output_struct, 'ranges') ...
                && isfield(output_struct, 'target_positions') ...
                && ~obj.visible
                obj.video_trajectory(...
                    output_struct.scatterer_positions, ...
                    output_struct.t_m, ...
                    output_struct.ranges, ...
                    output_struct.target_positions)
            end

            if (obj.bool_plot_range || obj.bool_plot_all) ...
                    && isfield(output_struct, 'ranges') ...
                    && isfield(output_struct, 't_m')
                obj.plot_range(output_struct.ranges, ...
                    output_struct.t_m)
            end

            if (obj.bool_plot_raw_isar || obj.bool_plot_all) ...
                    && isfield(output_struct, 'rx_signal') ...
                    && isfield(output_struct, 'range_array') ...
                    && isfield(output_struct, 'ranges')
                obj.plot_raw_isar(output_struct.rx_signal, ...
                    output_struct.range_array, ...
                    output_struct.ranges)
            end

            if (obj.bool_plot_range_compressed || obj.bool_plot_all) ...
                    && isfield(output_struct, 'rx_signal_range_compressed')
                obj.plot_range_compressed(output_struct.rx_signal_range_compressed)
            end

            if (obj.bool_plot_range_tracking || obj.bool_plot_all) ...
                    && isfield(output_struct, 'rx_signal_aligned')
                obj.plot_range_tracking(output_struct.rx_signal_aligned)
            end

            if (obj.bool_plot_autofocus || obj.bool_plot_all) ...
                    && isfield(output_struct, 'rx_autofocused')
                obj.plot_autofocus(output_struct.rx_autofocused)
            end

            if (obj.bool_plot_rd || obj.bool_plot_all) ...
                    && isfield(output_struct, 'rx_signal_rd')
                obj.plot_rd(output_struct.rx_signal_rd)
            end

            if (obj.bool_plot_target || obj.bool_plot_all) ...
                    && isfield(output_struct, 'ranges') ...
                    && isfield(output_struct, 'target')
                obj.plot_target(output_struct.ranges, ...
                    output_struct.target)
            end

        end

        function plot_rd(obj, rx_signal_rd)

            % find point target
            [~, max_idx] = max(abs(rx_signal_rd), [], "all");
            [max_row, max_col] = ind2sub(size(rx_signal_rd), max_idx);
            
            if obj.visible
                figure
            else
                f = figure('Visible','off');
            end
            subplot(1,2,1)
            imagesc(obj.range_array, 1:size(rx_signal_rd, 1), abs(rx_signal_rd))
            title('Range-Doppler ISAR image', 'FontSize', obj.title_font_size)
            xlabel('Range [m]', 'FontSize', obj.label_font_size)
            ylabel('Doppler bins', 'FontSize',obj.label_font_size)
            colorbar
            axis square
            
            subplot(1,2,2)
            semilogy(abs(rx_signal_rd(:, max_col)))
            title('Doppler frequency for range with dominant scatterer', 'FontSize', obj.title_font_size)
            ylabel('Log scale amplitude', 'FontSize', obj.label_font_size)
            xlabel('Doppler bins', 'FontSize', obj.label_font_size)
            grid on
            axis square
            
            buffers = [25,10,5,0];
            for ibuffer = 1:size(buffers, 2)
                
                buffer = buffers(ibuffer);
                row_idx = max_row-buffer:max_row+buffer;
                col_idx = max_col-buffer:max_col+buffer;
                
                if col_idx(end) < size(obj.range_array,1)
                    break
                end
            end

            if obj.visible
                figure
            else
                f = figure('Visible','off');
            end
            imagesc(obj.range_array(col_idx), row_idx, ...
                abs(rx_signal_rd(row_idx, col_idx)))
            title('Range-Doppler ISAR image - Zoomed near dominant scatterer','FontSize', obj.title_font_size)
            xlabel('Range [m]', 'FontSize', obj.label_font_size)
            ylabel('Doppler bins', 'FontSize', obj.label_font_size)
            colorbar
            axis square

            if ~obj.visible
                set(gcf, 'Position', get(0, 'Screensize'));
                saveas(f, 'plots/rd.png')
            end

        end

        function plot_autofocus(obj, rx_autofocused)

            if obj.visible
                figure
            else
                f = figure('Visible','off');
            end
            % subplot(1,2,1)
            imagesc(obj.range_array, 1:size(rx_autofocused, 1), abs(rx_autofocused))
            title('Absolute value of ISAR image post autofocus', 'FontSize', obj.title_font_size)
            xlabel('range [m]', 'FontSize', obj.label_font_size)
            ylabel('pulse index', 'FontSize', obj.label_font_size)
            colorbar
            axis square
            
            % subplot(1,2,2)
            % imagesc(obj.range_array, 1:size(rx_autofocused, 1), real(rx_autofocused))
            % title('Real of ISAR image post autofocus')
            % xlabel('range [m]')
            % ylabel('pulse index')
            % colorbar
            % axis square

            if ~obj.visible
                set(gcf, 'Position', get(0, 'Screensize'));
                saveas(f, 'plots/autofocus.png')
            end

        end

        function plot_range_tracking(obj, rx_signal_aligned)

            if obj.visible
                figure
            else
                f = figure('Visible','off');
            end
            % subplot(1,2,1)
            imagesc(obj.range_array, 1:size(rx_signal_aligned, 1), abs(rx_signal_aligned))
            title('Absolute value of ISAR image post range tracking', 'FontSize', obj.title_font_size)
            xlabel('range [m]', 'FontSize', obj.label_font_size)
            ylabel('pulse index', 'FontSize', obj.label_font_size)
            colorbar
            axis square
            
            % subplot(1,2,2)
            % imagesc(obj.range_array, 1:size(rx_signal_aligned, 1),real(rx_signal_aligned))
            % title('Real of ISAR image post range tracking')
            % xlabel('range [m]')
            % ylabel('pulse index')
            % colorbar
            % axis square

            if ~obj.visible
                set(gcf, 'Position', get(0, 'Screensize'));
                saveas(f, 'plots/range_tracking.png')
            end

        end

        function plot_range_compressed(obj, rx_signal_range_compressed)

            if obj.visible
                figure
            else
                f = figure('Visible','off');
            end
            % subplot(1,2,1)
            imagesc(obj.range_array, 1:size(rx_signal_range_compressed, 1), abs(rx_signal_range_compressed))
            title('Absolute value range compressed ISAR image', 'FontSize', obj.title_font_size)
            xlabel('range [m]', 'FontSize', obj.label_font_size)
            ylabel('pulse index', 'FontSize', obj.label_font_size)
            colorbar
            axis square
            
            % subplot(1,2,2)
            % imagesc(obj.range_array, 1:size(rx_signal_range_compressed, 1), real(rx_signal_range_compressed))
            % title('Real Range compressed ISAR image')
            % xlabel('range [m]')
            % ylabel('pulse index')
            % colorbar
            % axis square

            if ~obj.visible
                set(gcf, 'Position', get(0, 'Screensize'));
                saveas(f, 'plots/range_compression.png')
            end

        end

        function plot_raw_isar(obj, rx_signal, range_array, ranges)

            if obj.visible
                figure
            else
                f = figure('Visible','off');
            end
            % subplot(1,2,1)
            % imagesc(obj.range_array, 1:size(rx_signal, 1), imag(rx_signal))
            % title('Imaginary value of received signal - raw ISAR data')
            % xlabel('range [m]')
            % ylabel('pulse index')
            % colorbar
            % axis square
            
            % subplot(1,2,1)
            imagesc(obj.range_array, 1:size(rx_signal, 1), real(rx_signal))
            title('Real value of received signal - raw ISAR data', 'FontSize', obj.title_font_size)
            xlabel('range [m]', 'FontSize', obj.label_font_size)
            ylabel('pulse index', 'FontSize', obj.label_font_size)
            colorbar
            axis square
            
            if ~obj.visible
                set(gcf, 'Position', get(0, 'Screensize'));
                saveas(f, 'plots/raw_data1.png')
            end
            % 
            if obj.visible
                figure
            else
                f = figure('Visible','off');
            end
            % subplot(1,2,1)
            % plot(range_array, imag(rx_signal(round(size(ranges,1)/2), :)))
            % hold on
            % xline(ranges(round(size(ranges,1)/2)), 'Color', 'r', 'LineWidth', 2, 'LineStyle','--')
            % title('Imaginary value of received signal for a singular pulse - raw ISAR data')
            % xlabel('range [m]')
            % ylabel('Amplitude')

            % subplot(1,2,2)
            plot(range_array, real(rx_signal(round(size(ranges,1)/2), :)))
            hold on
            xline(ranges(round(size(ranges,1)/2)), 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--')
            title('Real value of received signal for a singular pulse - raw ISAR data', 'FontSize', obj.title_font_size)
            xlabel('range [m]', 'FontSize', obj.label_font_size)
            ylabel('Amplitude', 'FontSize', obj.label_font_size)
            axis square


            if ~obj.visible
                set(gcf, 'Position', get(0, 'Screensize'));
                saveas(f, 'plots/raw_data2.png')
            end
            
            if obj.visible
                figure
            else
                f = figure('Visible','off');
            end
            s = surf(real(rx_signal));
            s.EdgeColor = 'none';
            xlabel('range bins')
            ylabel('pulses')

            if ~obj.visible
                set(gcf, 'Position', get(0, 'Screensize'));
                saveas(f, 'plots/raw_data3.png')
            end

        end

        function plot_range(obj, ranges, t_m)

            if obj.visible
                figure
            else
                f = figure('Visible','off');
            end
            plot(t_m, ranges)
            title('Target Range')
            xlabel('Slow time')
            ylabel('Range')
            axis square

            if ~obj.visible
                set(gcf, 'Position', get(0, 'Screensize'));
                saveas(f, 'plots/range.png')
            end

        end

        function plot_target(obj, ranges, target)

            if obj.visible
                figure
            else
                f = figure('Visible','off');
            end

            target_config = ...
                target.scatter_configuration(:,1:2) ...
                + ranges(1);

            plot(target_config(:,1), target_config(:,2), 'Marker','*', 'LineStyle','none', 'MarkerSize', 10)
            title('Target Scatter Initial positions')
            xlabel('x (m)')
            ylabel('y (m)')

            buffer = 10;
            xlim([min(target_config(:,1))-buffer, max(target_config(:,1))+buffer])
            ylim([min(target_config(:,2))-buffer, max(target_config(:,2))+buffer])
            grid on
            axis square

            if ~obj.visible
                set(gcf, 'Position', get(0, 'Screensize'));
                saveas(f, 'plots/target.png')
            end

        end

        function video_trajectory(obj, scatterer_positions, t_m, ranges, target_positions)

            % target positions dimensions:
            % [nRanges x nScatters x Positions]

            f = figure('Visible','off');
            set(gcf, 'Position', get(0, 'Screensize'));
            subplot(1,2,1)
            xmin = min([scatterer_positions(:,:,1); zeros(1, size(scatterer_positions, 2))], [],"all");
            ymin = min([scatterer_positions(:,:,2); zeros(1, size(scatterer_positions, 2))], [], "all");
            zmin = min([scatterer_positions(:,:,3); zeros(1, size(scatterer_positions, 2))], [], "all");

            
            plot3( ...
                target_positions(:,1), ...
                target_positions(:,2), ...
                target_positions(:,3), 'DisplayName', 'Target Trajectory', ...
                'LineWidth', 2)
            
            hold on

            xlabel('X [m]')
            ylabel('Y [m]')
            zlabel('Z [m]')
            xlim([xmin, max(scatterer_positions(:,:,1:2),[], "all")])
            ylim([ymin, max(scatterer_positions(:,:,1:2),[], "all")])
            zlim([zmin, max(scatterer_positions(:,:,3), [], "all")])
            view([0 45])
            grid on

            subplot(1,2,2)
            plot(t_m, ranges, 'LineWidth', 2)
            title('Target Range')
            xlabel('Slow time [s]')
            ylabel('Range [m]')
            axis square
            grid on
            hold on

            % we want a maximum of 500 images, if the ranges
            % exceed 500 - downsample
            total_frames = 50;
            if size(scatterer_positions,1) > total_frames
                step = ...
                    round(size(scatterer_positions,1) / total_frames);
                idx = 1:step:size(scatterer_positions,1);
                step10 = round(size(idx,2) / 10);
                idx10 = 1:step10:size(idx,2);
                idx10 = idx(idx10);
            else
                idx = 1:size(scatterer_positions,1);
            end
        
            j = 1;
            fprintf('Creating video of target trajectory.\n')
            for irange = idx

                if ismember(irange, idx10)
                    fprintf('   %3.1f percent complete\n', (irange*100/size(scatterer_positions,1)))
                end

                los_x = linspace(0, target_positions(irange, 1), 20);
                los_y = linspace(0, target_positions(irange, 2), 20);
                los_z = linspace(0, target_positions(irange, 3), 20);
                
                subplot(1,2,1)
                h = plot3(los_x, los_y, los_z, ...
                    'Color', 'r', 'DisplayName', 'LOS to target center', 'LineWidth', 2, 'LineStyle','--');
                hold on
                g = [];
                for iscatter = 1:size(scatterer_positions, 2)
                
                    
                    g = [g, plot3(scatterer_positions(irange, iscatter, 1), ...
                        scatterer_positions(irange, iscatter, 2), ...
                        scatterer_positions(irange, iscatter, 3), ...
                    'Color', 'r', 'DisplayName', 'Scatterers', 'Marker','*', 'MarkerSize', 5, 'LineStyle','none')];
                
                    
                end

                ff = fill3(scatterer_positions(irange, :, 1), ...
                    scatterer_positions(irange, :, 2), ...
                    scatterer_positions(irange, :, 3), 'r');

                legend([h, g(1)], {'LOS', 'scatterers'},'Location','northeast')
                k = [];
                for iscatter = 1:size(scatterer_positions, 2)
                    subplot(1,2,2)
                    k = [k, plot(t_m(irange), ...
                        ranges(irange, iscatter), '*', 'Color', 'r', 'MarkerSize', 10, 'Marker','o')];
                end

                F(j) = getframe(gcf);
                j = j + 1;
                delete(h)
                delete(k)
                delete(g)
                delete(ff)
            end

            % create the video writer with 1 fps
            writerObj = VideoWriter('plots/target_trajectory', 'MPEG-4');

            writerObj.FrameRate = 10; 

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

        end

        function plot_trajectory(obj, scatterer_positions)
            
            xmin = min([scatterer_positions(:,1); 0]);
            ymin = min([scatterer_positions(:,2); 0]);
            zmin = min([scatterer_positions(:,3); 0]);

            if obj.visible
                figure
            else
                f = figure('Visible','off');
            end
            
            plot3( ...
                scatterer_positions(:,1), ...
                scatterer_positions(:,2), ...
                scatterer_positions(:,3), 'DisplayName', 'Target Trajectory')
            
            hold on

            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            xlim([xmin, max(scatterer_positions(:,1:2),[], "all")])
            ylim([ymin, max(scatterer_positions(:,1:2),[], "all")])
            zlim([zmin, max(scatterer_positions(:,3))])
            
            ipt = 1;

            los_x = linspace(0, scatterer_positions(ipt,1), 20);
            los_y = linspace(0, scatterer_positions(ipt,2), 20);
            los_z = linspace(0, scatterer_positions(ipt,3), 20);
            h = plot3(los_x, los_y, los_z, 'Color', 'r', 'DisplayName', 'LOS');
            legend('Location','northeast')

            if ~obj.visible
                set(gcf, 'Position', get(0, 'Screensize'));
                saveas(f, 'plots/trajectory.png')
            end
                
        end
    end
end

