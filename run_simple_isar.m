
%% define the scenario
 
% instantiate constants
const = Constants;
c = const.c;

% simulation duration
T = 1; % [s]
            
% initialize radar parameters and LFM signal model
fc  = 10 * const.GHz2Hz; % [Hz] center frequency - X-band
B   = 149.9 * const.MHz2Hz; % [Hz] bandwidth
prf = 1 * const.kHz2Hz; % [Hz] pulse repetition frequency
fs  = 300 * const.MHz2Hz; % [Hz] sampling frequency
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

t_tx = (-Tp/2 : dt_fast_time : Tp/2 - dt_fast_time)';
tx_signal = sinc(B * t_tx); % baseband


%% define the target scatterers

% number of scatterers
num_scatterers = 3;

% define the shape of the simulated aircraft
% as a triangle relative to the center
radius = [3, 7];

% define the target's initial position
target_center_position = [0, 1000, 0];

scatter_relative_positions = ...
    [0, 0, 0; ... % scatter 1
    radius(2), 0, 0; ... % scatter 2
    0, radius(1), 0];     % scatter 3

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
fprintf('Propagating simulation\n')

% iterate through each pulse (column)
for ipulse = 1:num_pulses

    if mod(ipulse, num_pulses/10) == 0
    fprintf('   %i percent complete\n', (ipulse*100/num_pulses))
    end
    
    % iterate through each target scatter
    for ipt = 1:num_scatterers
    
        % extract the position of each scatterer
        % relative to the axis of rotation
        x = scatter_relative_positions(ipt,1); % [m]
        y = scatter_relative_positions(ipt,2); % [m]
        R0 = norm(target_center_position); % [m]
        theta = yaws(ipulse);
        
        R = [cos(yaw), -sin(yaw), 0; ...
                 sin(yaw), cos(yaw), 0; ...
                 0, 0, 1];
        
        scatterer_absolute_position = ...
            (R * scatter_relative_positions(ipt,:)')' ...
            + target_center_position;
        range1 = norm(scatterer_absolute_position);

        % calculate the range of the scatterer relative
        % to the radar at the origin (0,0)
        range = R0 + x*sin(theta) + y*cos(theta); % [m]
        delay = 2 * range / c; % [s]
        
        % define the reflectivity of the scatterer
        A = 1;
        
        % define the phase 
        phase = exp( -1j * 4 * pi * fc * range / c);
        
        % define the echo, equation 2 - not sure if this
        % should be normalized or unnormalized
        scatterer_signal = ...
            A * sinc(B * (t_fast - delay)) * phase;
        
        % add the scatterer's contribution to the range
        % profile
        rx_signal(ipulse, :) = ...
            rx_signal(ipulse, :) + scatterer_signal';
        
        % save target position data
        ranges(ipulse,ipt) = range;
        R = [cos(yaw), -sin(yaw), 0; ...
                 sin(yaw), cos(yaw), 0; ...
                 0, 0, 1];

        scatterer_position = ...
            R * scatter_relative_positions(ipt,:)';
        scatterer_positions(ipulse, ipt, :) = ...
            scatterer_position';
    end
end


%% range compression

% range compression via match filtering using the
% transmitted signal pulse
h = conj(flipud(tx_signal));
rx_signal_range_compressed = zeros(size(rx_signal));
for ipulse =1:num_pulses    
    rx_signal_range_compressed(ipulse,:) = ...
        conv(rx_signal(ipulse,:), h, "same");
end

% adjust the range axis by the match filter's group delay
% range_array = range_array - (Tp * c / 2);


%% backprojection

% define the cross range resolution for the grid
% cross_range_resolution =  c / (2 * yaws(end) * fc);
cross_range_resolution = 0.5;

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
    
% f = waitbar(0, 'Performing backprojection');
% iterate through each pixel
fprintf('Performing backprojection\n')
for ix = 1:Nx % cross range
    
    if mod(ix, round(Nx/10)) == 0
        % waitbar(ix/round(Nx),f,sprintf('Progress: %d %%', floor(ix/Nx*100)))
        fprintf('   %3.1f percent complete: %4.2f seconds\n',...
            (ix*100/round(Nx)), toc(t))
        t = tic;
    end

    for iy = 1:Ny % range

        % extract the pixel location relative to the
        % center of the target
        pixel_location = [x_array(ix), y_array(iy), 0];

        % initialize the image value as zero
        image_value = 0;

        % add a hamming window to surpress sidelobes
        window = hamming(num_pulses);

        % iterate through each pulse
        for ipulse = 1:num_pulses

            % extract yaw
            yaw = yaws(ipulse);

            % calculate the location of the pixel
            % relative to the radar
            R = [cos(yaw), -sin(yaw), 0; ...
                 sin(yaw), cos(yaw), 0; ...
                 0, 0, 1];

            pixel_location_radar = ...
                (R * pixel_location')' ...
                + target_center_position;

            % find the range of the pixel to the radar
            range = norm(pixel_location_radar);

            % interpolate the range
            echo = interp1(...
                range_array, ...
                rx_signal_range_compressed(ipulse, :), ...
                range, ...
                'pchip', 0);

            % Apply phase correction (match propagation delay)
            phase = ...
                exp(1j * 4 * pi * fc * range / c);

            % calculation the pulse contribution
            pulse_contribution = ...
                window(ipulse) * echo * phase;

            % Sum contribution
            image_value = ...
                image_value + pulse_contribution;

        end
        rx_signal_bp(ix,iy) = image_value; % [cross range x range]

    end
end


%% trajectory

% range and trajectory

figure
subplot(1,2,1)
plot(t_slow, ranges, 'LineWidth', 2)
title('Target Range', 'FontSize', 24)
xlabel('Slow time', 'FontSize', 16)
ylabel('Range', 'FontSize', 16)
axis square
ax = gca;
set(ax,'FontSize',16)
ax.YDir = "reverse";

xmin = min([scatterer_positions(:,:,1)], [], "all");
ymin = min([scatterer_positions(:,:,2)], [], "all");
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

title('Target Trajectory', 'FontSize', 24)
xlabel('X [m]', 'FontSize', 16)
ylabel('Y [m]', 'FontSize', 16)
grid on
axis square
ax = gca;
set(ax,'FontSize',16)
ax.YDir = "reverse";

% range compression

% figure
% subplot(1,2,1)
% imagesc(range_array, 1:size(rx_signal_range_compressed, 1), abs(rx_signal_range_compressed))
% title('Absolute value range compressed ISAR image', 'FontSize', 24)
% xlabel('range [m]', 'FontSize', 16)
% ylabel('pulse index', 'FontSize', 16)
% colorbar
% axis square
% set(gca,'FontSize',16)
% 
% pulse_idx = round(size(rx_signal_range_compressed,1)/2);
% subplot(1,2,2)
% plot(range_array, log(abs(rx_signal_range_compressed(pulse_idx, :))))
% title(['Range profile for a singular range compressed pulse (' num2str(pulse_idx), ')'], 'FontSize', 24)
% xlabel('range [m]', 'FontSize', 16)
% ylabel('Log of Signal', 'FontSize', 16)
% colorbar
% axis square
% set(gca,'FontSize',16)

% backprojection

figure
subplot(1,2,1)
imagesc(x_array, y_array, 20*log10(abs(rx_signal_bp.') + eps));
title('Backprojection image - Log scaled', 'FontSize', 24)
xlabel('x (cross-range)', 'FontSize', 16)
ylabel('y (range)', 'FontSize', 16)
axis square
colorbar   
set(gca,'FontSize',16)

subplot(1,2,2)
imagesc(x_array, y_array, abs(rx_signal_bp.'))
title('Backprojection image', 'FontSize', 24)
xlabel('x (cross-range)', 'FontSize', 16)
ylabel('y (range)', 'FontSize', 16)
axis square
colorbar   
set(gca,'FontSize',16)
            







