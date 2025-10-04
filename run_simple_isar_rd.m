
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
doppler_array = linspace(-prf/2, prf/2, size(t_slow,1))';

% calculate the range resolution
range_resolution = c / (2 * B);

% calculate the number of pulses
num_pulses = round(T / dt_slow);
num_range_bins = size(t_fast,1);


%% define the target scatterers

% number of scatterers
num_scatterers = 1;

% define the target's initial position
target_center_position = [0, 1000, 0];

scatter_relative_positions = [7, 0, 0];

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
            rx_signal(ipulse, :) + scatterer_signal';
        
        % save target position data
        ranges(ipulse,ipt) = range;

        scatterer_position = ...
            R * scatter_relative_positions(ipt,:)';
        scatterer_positions(ipulse, ipt, :) = ...
            scatterer_position';
    end
end


%% rd processing

% perform Doppler processing (or az FFT) on the range
% compressed data   
rx_signal_rd = ...
    fftshift(fft(rx_signal, [], 1), 1);


%% plotting

% range and trajectory

f = figure('Visible','on');
subplot(1,2,1)
plot(t_slow, ranges, 'LineWidth', 2)
title('Target Range', 'FontSize', 24)
xlabel('Slow time', 'FontSize', 16)
ylabel('Range', 'FontSize', 16)
axis square
grid on
ax = gca;
set(ax,'FontSize',16)
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

title('Target Trajectory', 'FontSize', 24)
xlabel('X [m]', 'FontSize', 16)
ylabel('Y [m]', 'FontSize', 16)
grid on
axis square
ax = gca;
set(ax,'FontSize',16)
ax.YDir = "reverse";
set(gcf, 'Position', get(0, 'Screensize'));
%saveas(f,'plots/range_and_trj.png')

% range compression

f = figure('Visible','on');
subplot(1,2,1)
imagesc(range_array, 1:size(rx_signal, 1), abs(rx_signal))
title('Absolute value range compressed ISAR image', 'FontSize', 24)
xlabel('range [m]', 'FontSize', 16)
ylabel('pulse index', 'FontSize', 16)
colorbar
axis square
set(gca,'FontSize',16)
xlim([990,1010])

pulse_idx = round(size(rx_signal,1)/2);
subplot(1,2,2)
plot(range_array, log(abs(rx_signal(pulse_idx, :))))
title(['Range profile for a singular range compressed pulse (' num2str(pulse_idx), ')'], 'FontSize', 24)
xlabel('range [m]', 'FontSize', 16)
ylabel('log-scaled signal', 'FontSize', 16)
colorbar
axis square
set(gca,'FontSize',16)
set(gcf, 'Position', get(0, 'Screensize'));
xlim([990,1010])
%saveas(f,'plots/range_compression.png')


% range-doppler

f = figure('Visible','on');
subplot(1,2,1)
imagesc(doppler_array, range_array, 20*log10(abs(rx_signal_rd.')));
title('Range-Doppler Processing image - Log scaled', 'FontSize', 24)
ylabel('Range [m]', 'FontSize', 16)
xlabel('Doppler [Hz]', 'FontSize',16)
axis square
grid on
colorbar   
set(gca,'FontSize',16)

subplot(1,2,2)
imagesc(doppler_array, range_array, abs(rx_signal_rd.'));
title('Range-Doppler Processing', 'FontSize', 24)
ylabel('Range [m]', 'FontSize', 16)
xlabel('Doppler [Hz]', 'FontSize',16)
axis square
grid on
colorbar   
set(gca,'FontSize',16)

%saveas(f, 'plots/rd.png')






