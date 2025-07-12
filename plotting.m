classdef plotting < handle
    %PLOTTING Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        plot_range_bool = false
        plot_target_trajectory_bool = false
        video_target_trajectory_bool = false

    end
    
    methods

        function obj = plotting()
            
        end

        function plot(obj, output_struct)
            if obj.plot_target_trajectory_bool
                obj.plot_trajectory(output_struct.target_positions)
            end
            if obj.video_target_trajectory_bool
                obj.video_trajectory(output_struct.target_positions)
            end
        end

        function plot_range(obj, ranges)

            figure
            plot(ranges)
            title('Target Range')
            xlabel('Slow time')
            ylabel('Range')

        end

        function video_trajectory(obj, target_positions)

            f = figure('Visible','off');
            
            plot3( ...
                target_positions(:,1), ...
                target_positions(:,2), ...
                target_positions(:,3), 'DisplayName', 'Target Trajectory')
            
            hold on

            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            xlim([0, max(target_positions(:,1:2),[], "all")])
            ylim([0, max(target_positions(:,1:2),[], "all")])
            zlim([0, max(target_positions(:,3))])
            
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

            writerObj.FrameRate = 1; % set the seconds per image

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

            f = figure('Visible','on');
            
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
                
        end
    end
end

