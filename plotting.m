classdef plotting < handle
    %PLOTTING Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        visible = false
        plot_all = false

        plot_range_bool = false
        plot_target_trajectory_bool = false
        video_target_trajectory_bool = false

        plot_raw_isar_bool = false
        plot_range_compressed_bool = false
        plot_range_tracking_bool = false
        plot_autofocus_bool = false
        plot_rd_bool = false

        range_array

    end
    
    methods

        function obj = plotting()
            
        end

        function plot(obj, output_struct)

            if obj.plot_range_bool ...
                || obj.plot_target_trajectory_bool ...
                || obj.video_target_trajectory_bool ...
                || obj.plot_raw_isar_bool ...
                || obj.plot_range_compressed_bool ...
                || obj.plot_range_tracking_bool ...
                || obj.plot_autofocus_bool ...
                || obj.plot_rd_bool

                delete plots/*.png
            end

            if isfield(output_struct, 'range_array')
                obj.range_array = output_struct.range_array;
            end

            if (obj.plot_target_trajectory_bool || obj.plot_all) ...
                    && isfield(output_struct, 'target_positions')
                obj.plot_trajectory(output_struct.target_positions)
            end

            if (obj.video_target_trajectory_bool || obj.plot_all)...
                && isfield(output_struct, 'target_positions') ...
                && ~obj.visible
                obj.video_trajectory(output_struct.target_positions)
            end

            if (obj.plot_range_bool || obj.plot_all) ...
                    && isfield(output_struct, 'ranges') ...
                    && isfield(output_struct, 't_m')
                obj.plot_range(output_struct.ranges, ...
                    output_struct.t_m)
            end

            if (obj.plot_raw_isar_bool || obj.plot_all) ...
                    && isfield(output_struct, 'rx_signal') ...
                    && isfield(output_struct, 'range_array') ...
                    && isfield(output_struct, 'ranges')
                obj.plot_raw_isar(output_struct.rx_signal, ...
                    output_struct.range_array, ...
                    output_struct.ranges)
            end

            if (obj.plot_range_compressed_bool || obj.plot_all) ...
                    && isfield(output_struct, 'rx_signal_range_compressed')
                obj.plot_range_compressed(output_struct.rx_signal_range_compressed)
            end

            if (obj.plot_range_tracking_bool || obj.plot_all) ...
                    && isfield(output_struct, 'rx_signal_aligned')
                obj.plot_range_tracking(output_struct.rx_signal_aligned)
            end

            if (obj.plot_autofocus_bool || obj.plot_all) ...
                    && isfield(output_struct, 'rx_autofocused')
                obj.plot_autofocus(output_struct.rx_autofocused)
            end

            if (obj.plot_rd_bool || obj.plot_all) ...
                    && isfield(output_struct, 'rx_signal_rd')
                obj.plot_rd(output_struct.rx_signal_rd)
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
            title('Range-Doppler ISAR image')
            xlabel('Range [m]')
            ylabel('Doppler bins')
            colorbar
            axis square
            
            subplot(1,2,2)
            semilogy(abs(rx_signal_rd(:, max_col)))
            title('Doppler frequency for range with dominant scatterer')
            ylabel('Log scale amplitude')
            xlabel('Doppler bins')
            grid on
            axis square
            
            buffer = 25;
            
            row_idx = max_row-buffer:max_row+buffer;
            col_idx = max_col-buffer:max_col+buffer;
            
            if obj.visible
                figure
            else
                f = figure('Visible','off');
            end
            imagesc(obj.range_array(col_idx), row_idx, ...
                abs(rx_signal_rd(row_idx, col_idx)))
            title('Range-Doppler ISAR image - Zoomed near dominant scatterer')
            xlabel('Range [m]')
            ylabel('Doppler bins')
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
            subplot(1,2,1)
            imagesc(obj.range_array, 1:size(rx_autofocused, 1), abs(rx_autofocused))
            title('Absolute value of ISAR image post autofocus')
            xlabel('range [m]')
            ylabel('pulse index')
            colorbar
            axis square
            
            subplot(1,2,2)
            imagesc(obj.range_array, 1:size(rx_autofocused, 1), real(rx_autofocused))
            title('Real of ISAR image post autofocus')
            xlabel('range [m]')
            ylabel('pulse index')
            colorbar
            axis square

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
            subplot(1,2,1)
            imagesc(obj.range_array, 1:size(rx_signal_aligned, 1), abs(rx_signal_aligned))
            title('Absolute value of ISAR image post range tracking')
            xlabel('range [m]')
            ylabel('pulse index')
            colorbar
            axis square
            
            subplot(1,2,2)
            imagesc(obj.range_array, 1:size(rx_signal_aligned, 1),real(rx_signal_aligned))
            title('Real of ISAR image post range tracking')
            xlabel('range [m]')
            ylabel('pulse index')
            colorbar
            axis square

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
            subplot(1,2,1)
            imagesc(obj.range_array, 1:size(rx_signal_range_compressed, 1), abs(rx_signal_range_compressed))
            title('Absolute value range compressed ISAR image')
            xlabel('range [m]')
            ylabel('pulse index')
            colorbar
            axis square
            
            subplot(1,2,2)
            imagesc(obj.range_array, 1:size(rx_signal_range_compressed, 1), real(rx_signal_range_compressed))
            title('Real Range compressed ISAR image')
            xlabel('range [m]')
            ylabel('pulse index')
            colorbar
            axis square

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
            subplot(1,2,1)
            imagesc(obj.range_array, 1:size(rx_signal, 1), imag(rx_signal))
            title('Imaginary value of received signal - raw ISAR data')
            xlabel('range [m]')
            ylabel('pulse index')
            colorbar
            
            subplot(1,2,2)
            imagesc(obj.range_array, 1:size(rx_signal, 1), real(rx_signal))
            title('Real value of received signal - raw ISAR data')
            xlabel('range [m]')
            ylabel('pulse index')
            colorbar
            
            if ~obj.visible
                set(gcf, 'Position', get(0, 'Screensize'));
                saveas(f, 'plots/raw_data1.png')
            end

            if obj.visible
                figure
            else
                f = figure('Visible','off');
            end
            subplot(1,2,1)
            plot(range_array, imag(rx_signal(round(size(ranges,1)/2), :)))
            hold on
            xline(ranges(round(size(ranges,1)/2)), 'Color', 'r', 'LineWidth', 2, 'LineStyle','--')
            title('Imaginary value of received signal for a singular pulse - raw ISAR data')
            xlabel('range [m]')
            ylabel('Amplitude')

            subplot(1,2,2)
            plot(range_array, real(rx_signal(round(size(ranges,1)/2), :)))
            hold on
            xline(ranges(round(size(ranges,1)/2)), 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--')
            title('Real value of received signal for a singular pulse - raw ISAR data')
            xlabel('range [m]')
            ylabel('Amplitude')


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

            if ~obj.visible
                set(gcf, 'Position', get(0, 'Screensize'));
                saveas(f, 'plots/range.png')
            end

        end

        function video_trajectory(obj, target_positions)

            f = figure('Visible','off');
            
            xmin = min([target_positions(:,1); 0]);
            ymin = min([target_positions(:,2); 0]);
            zmin = min([target_positions(:,3); 0]);

            plot3( ...
                target_positions(:,1), ...
                target_positions(:,2), ...
                target_positions(:,3), 'DisplayName', 'Target Trajectory')
            
            hold on

            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            xlim([xmin, max(target_positions(:,1:2),[], "all")])
            ylim([ymin, max(target_positions(:,1:2),[], "all")])
            zlim([zmin, max(target_positions(:,3))])
            
            for ipt = 1:size(target_positions(:,1))

                los_x = linspace(0, target_positions(ipt,1), 20);
                los_y = linspace(0, target_positions(ipt,2), 20);
                los_z = linspace(0, target_positions(ipt,3), 20);
                h = plot3(los_x, los_y, los_z, 'Color', 'r', 'DisplayName', 'LOS');
                legend('Location','northeast')
                F(ipt) = getframe(gcf);
                delete(h)
            end

            % create the video writer with 1 fps
            writerObj = VideoWriter('plots/target_trajectory', 'MPEG-4');

            writerObj.FrameRate = 60; 

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

        function plot_trajectory(obj, target_positions)
            
            xmin = min([target_positions(:,1); 0]);
            ymin = min([target_positions(:,2); 0]);
            zmin = min([target_positions(:,3); 0]);

            if obj.visible
                figure
            else
                f = figure('Visible','off');
            end
            
            plot3( ...
                target_positions(:,1), ...
                target_positions(:,2), ...
                target_positions(:,3), 'DisplayName', 'Target Trajectory')
            
            hold on

            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            xlim([xmin, max(target_positions(:,1:2),[], "all")])
            ylim([ymin, max(target_positions(:,1:2),[], "all")])
            zlim([zmin, max(target_positions(:,3))])
            
            ipt = 1;

            los_x = linspace(0, target_positions(ipt,1), 20);
            los_y = linspace(0, target_positions(ipt,2), 20);
            los_z = linspace(0, target_positions(ipt,3), 20);
            h = plot3(los_x, los_y, los_z, 'Color', 'r', 'DisplayName', 'LOS');
            legend('Location','northeast')

            if ~obj.visible
                set(gcf, 'Position', get(0, 'Screensize'));
                saveas(f, 'plots/trajectory.png')
            end
                
        end
    end
end

