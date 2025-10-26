
%% define the scenario

show_plots = true;
save_data = false;
save_plots = false;

% instantiate constants
const = Constants;
c = const.c;

% simulation duration
T = 1; % [s]
            
% initialize radar parameters and LFM signal model
fc  = 10 * const.GHz2Hz; % [Hz] center frequency - X-band
B   = 149.9 * const.MHz2Hz; % [Hz] bandwidth
prf = 1 * const.kHz2Hz; % [Hz] pulse repetition frequency
fs  = 600 * const.MHz2Hz; % [Hz] sampling frequency
Tp  = 10 * const.us2s; % [s] pulse width


%% basic signal model

% slow time, pulse response interval (PRI)
dt_slow = 1/prf; % [s]
t_slow = (0:dt_slow:T-dt_slow)';

% define the fast time array
dt_fast_time = 1/fs;
t_fast = (0:dt_fast_time:(Tp)-(1/fs))';

% define the chirp time array
range_array = t_fast .* c / 2;

% calculate the range resolution
range_resolution = c / (2 * B);

% calculate the number of pulses
num_pulses = round(T / dt_slow);
num_range_bins = size(t_fast,1);


%% define the target scatterers

% number of scatterers
num_scatterers = 3;

% define the target's initial position
target_center_position = [0, 1000, 0];

scatter_relative_positions = [7, 0, 0; ...
                              0, 3, 0; ...
                              0, 0, 0];

% calculate the yaws
yawing_rate = pi/4;
yaws = yawing_rate * t_slow;


%% raw isar image

% initialize received signal matrix
% [number of range cells x number of cross range cells]
rx_signal = zeros(num_pulses, num_range_bins);

% output vectors
scatterer_positions = zeros(num_pulses, num_scatterers, 3);
ranges = zeros(num_pulses, num_scatterers);
phases = zeros(num_pulses, num_scatterers);
fprintf('Propagating simulation\n')

% iterate through each pulse (column)
for ipulse = 1:num_pulses

    if mod(ipulse, num_pulses/10) == 0
    fprintf('   %i percent complete\n', (ipulse*100/num_pulses))
    end
    
    % iterate through each target scatter
    for ipt = 1:num_scatterers
    
        % formulate the rotation matrix about the z axis
        R = [cos(yaws(ipulse)), -sin(yaws(ipulse)), 0; ...
             sin(yaws(ipulse)), cos(yaws(ipulse)), 0; ...
             0, 0, 1];
        
        % calculate the scatterer's position relative to the
        % radar
        scatterer_absolute_position = ...
            (R * scatter_relative_positions(ipt,:)')' ...
            + target_center_position;
        range = norm(scatterer_absolute_position);

        % calculate the dealy
        delay = 2 * range / c; % [s]

        % define the reflectivity of the scatterer
        A = 1;
        
        % define the phase 
        phase = exp( -1j * 4 * pi * fc * range / c);
        
        % define the echo, equation 2 - not sure if this
        % should be normalized or unnormalized sinc
        scatterer_signal = ...
            A * sinc(B * (t_fast - delay)) * phase;
        
        % add the scatterer's contribution to the range
        % profile
        rx_signal(ipulse, :) = ...
            rx_signal(ipulse, :) + scatterer_signal.';
        
        % save target position data
        ranges(ipulse,ipt) = range;
        phases(ipulse,ipt) = phase;

        scatterer_position = ...
            R * scatter_relative_positions(ipt,:)';
        scatterer_positions(ipulse, ipt, :) = ...
            scatterer_position';
    end
end


%% backprojection

% define the cross range resolution for the grid
%cross_range_resolution =  c / (2 * yaws(end) * fc);
cross_range_resolution = 1;

% define the image dimensions using the range and
% cross-range resolution and the extent of the
% dimensions, 
cross_range_extent = 20; % [m]
range_extent = 20; % [m]

% the output of the backprojection algorithm will be 
% [Nx x Ny] where Nx is cross-range, Ny is range
Nx = round(cross_range_extent / cross_range_resolution);
Ny = round(range_extent / range_resolution);

% calculate the x and y values that define the image and
% the pixel locations relative to the radar
x_array = (-floor(Nx/2):ceil(Nx/2)-1) * cross_range_resolution;
y_array = (-floor(Ny/2):ceil(Ny/2)-1) * range_resolution;

% initialize final image
rx_signal_bp = zeros(Nx, Ny);

t = tic;

% formulate the mesh grid
[X,Y] = meshgrid(x_array,y_array);
[trow, tcol] = find(...
     scatter_relative_positions(1,1) == X ...
    & scatter_relative_positions(1,2) == Y);

% for debugging purposes calculate the difference between
% the range of target 2 and the pixel initially containing 
% target 2 (should be zero across every pulse if motion is 
% compensated for)
ranges_diff = zeros(num_pulses,1);
contributions = zeros(num_pulses, 1);
phases_bp = zeros(num_pulses, 1);

% iterate through each pixel
fprintf('Performing backprojection\n')
for ix = 1:Nx % cross range
    
    if mod(ix, round(Nx/10)) == 0
        fprintf('   %3.1f percent complete: %4.2f seconds\n',...
            (ix*100/round(Nx)), toc(t))
        t = tic;
    end

    for iy = 1:Ny % range

        % extract the pixel location relative to the
        % center of the target
        pixel_rel_location = [x_array(ix), y_array(iy), 0];

        % initialize the image value as zero
        image_value = 0;

        % iterate through each pulse
        for ipulse = 1:num_pulses

            % formulate the rotation matrix about the z axis
            R = [cos(yaws(ipulse)), -sin(yaws(ipulse)), 0; ...
                 sin(yaws(ipulse)), cos(yaws(ipulse)), 0; ...
                 0, 0, 1];

            pixel_location_radar = ...
                (R * pixel_rel_location')' ...
                + target_center_position;

            % find the range of the pixel to the radar
            range = norm(pixel_location_radar);

            % debugging
            if trow == iy && tcol == ix ...
                && ipulse == 500
                % find the range index
                [~,range_idx] = min(abs(range_array - range));

                % find the abs rx value corresponding to
                % this range
                abs_rx = abs(rx_signal(ipulse, range_idx));    

                % f = figure('Visible','on');
                
                % plot(range_array, 20*log10(abs(rx_signal(ipulse,:))), 'DisplayName', 'Log-scaled range profile', 'LineWidth', 2)
                % xline(range_array(range_idx), 'DisplayName', 'Range of rotated pixel', 'LineWidth', 1.5, 'Color', 'r')
                % yline(-13.3, 'DisplayName', 'Expected sidelobes (-13.3 dB)', 'Color', 'g', 'LineWidth', 1.5)
                % xlabel('Range [m]', 'FontSize', 16)
                % ylabel('Log-scaled range profile [dB]', 'FontSize', 16)
                % ax = gca;
                % set(ax,'FontSize',16)
                % legend
                % title('Range profile for a singular range compressed pulse (500)', 'FontSize', 24)
                % 
                % xlim([range_array(range_idx) - 20, range_array(range_idx) + 20])
                % set(gcf, 'Position', get(0, 'Screensize'));
                % saveas(f,'plots/backproj_range.png')

            end

            [~,minIdx] = min(abs(range - range_array));
            abs_rx = abs(rx_signal(ipulse,:));

            % interpolate the range
            echo = interp1(...
                range_array, ...
                rx_signal(ipulse, :), ...
                range, ...
                'linear', 0);

            % Apply phase correction (match propagation delay
            phase = exp( 1j * 4 * pi * fc * range / c);

            % calculation the pulse contribution
            pulse_contribution = echo * phase;

            % debugging
            if trow == iy && tcol == ix

                % save the range difference for the pixel
                % that should contain target 2
                ranges_diff(ipulse) = ...
                abs(range - ranges(ipulse,1));
                contributions(ipulse) = pulse_contribution;
                phases_bp(ipulse) = phase;
            end

            % Sum contribution
            image_value = image_value + pulse_contribution;

        end
        rx_signal_bp(ix,iy) = image_value; % [cross range x range]

    end
end


%% plotting
%close all
% total_contribution = sum(contributions); % complex sum
% ideal_magnitude = sum(abs(contributions)); % ideal if all in phase
% 
% efficiency = abs(total_contribution) / ideal_magnitude;
% fprintf('Phase alignment efficiency: %.2f%%\n', efficiency * 100);


% f = figure
% subplot(1,2,1)
% plot(angle(contributions), 'LineWidth', 2)
% xlabel('Pulses', 'FontSize', 16)
% ylabel('Phase [rad]', 'FontSize', 16)
% axis square
% grid on
% ax = gca;
% set(ax,'FontSize',16)
% subplot(1,2,2)
% plot(abs(contributions), 'LineWidth', 2)
% title('Magnitude at pixel where target starts at', 'FontSize', 24)
% xlabel('Pulses', 'FontSize', 16)
% ylabel('Magnitude', 'FontSize', 16)
% axis square
% grid on
% ax = gca;
% set(ax,'FontSize',16)
% title('Phase at pixel where target starts at', 'FontSize', 24)
% 
% set(gcf, 'Position', get(0, 'Screensize'));
% saveas(f,'plots/mag_and_phase.png')

% range and trajectory
visible = 'off';
if show_plots
    visible = 'on'; 
end

f = figure('Visible',visible);
subplot(1,2,1)
plot(t_slow, ranges, 'LineWidth', 2)
xlabel('Slow time', 'FontSize', 16)
ylabel('Range', 'FontSize', 16)
axis square
grid on
ax = gca;
set(ax,'FontSize',16)
title('Target Range', 'FontSize', 24)
ax.YDir = "reverse";

subplot(1,2,2)
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
title('Target Trajectory', 'FontSize', 24)

ax.YDir = "reverse";
set(gcf, 'Position', get(0, 'Screensize'));
if save_plots
    saveas(f,'plots/range_and_trj.png')
end
% range compression

f = figure('Visible',visible);
subplot(1,2,1)
title('Absolute value range compressed ISAR image', 'FontSize', 24)
imagesc(range_array, 1:size(rx_signal, 1), abs(rx_signal))
xlabel('range [m]', 'FontSize', 16)
ylabel('pulse index', 'FontSize', 16)
colorbar
axis square
set(gca,'FontSize',16)
xlim([990,1010])

pulse_idx = round(size(rx_signal,1)/2);
subplot(1,2,2)
plot(range_array, 20*log(abs(rx_signal(pulse_idx, :))))
title(['Range profile for a singular range compressed pulse (' num2str(pulse_idx), ')'], 'FontSize', 24)
xlabel('range [m]', 'FontSize', 16)
ylabel('log-scaled signal', 'FontSize', 16)
colorbar
axis square
set(gca,'FontSize',16)
set(gcf, 'Position', get(0, 'Screensize'));
xlim([990,1010])
if save_plots
    saveas(f,'plots/range_compression.png')
end

% backprojection

f = figure('Visible',visible);
subplot(1,2,1)
imagesc(x_array, y_array, 20*log10(abs(rx_signal_bp.') + eps));
xlabel('x (cross-range)', 'FontSize', 16)
ylabel('y (range)', 'FontSize', 16)
axis square
grid on
colorbar   
set(gca,'FontSize',16)
title('Backprojection image - Log scaled', 'FontSize', 24)


subplot(1,2,2)
imagesc(x_array, y_array, abs(rx_signal_bp.'))
xlabel('x (cross-range)', 'FontSize', 16)
ylabel('y (range)', 'FontSize', 16)
axis square
grid on
colorbar   
set(gca,'FontSize',16)
title('Backprojection image', 'FontSize', 24)

set(gcf, 'Position', get(0, 'Screensize'));
if save_plots
    saveas(f, 'plots/backproj.png')
end

% save data for compressed sensing use
if save_data
    save('data/rx_signal_bp_sd.mat', "rx_signal_bp");
end



