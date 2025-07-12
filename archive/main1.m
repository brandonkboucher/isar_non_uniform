
% The following is the summary I wrote after meeting with my
% advisors:

% Begin templating ISAR non-uniform sampling for an
% accelerating target resulting in quadratically spaced
% sampling. Simulate in MATLAB an observer on the ground
% and a target (begin with point) accelerating overhead
% with an angle of observation at each pulse. You can use
% an E&M simulator available like the one from Ansys/Nvidia.

% The following simulation follows the methodology described
% in: "Simulation of ISAR imaging of moving targets" by V.C.
% Chen and M.J. Miceli

%% Define the tx signal parameters

% instantiate constants
const = Constants;

% center frequency
fc = 16.7 * const.GHz2Hz; % [Hz]

% bandwidth
B = 150 * const.MHz2Hz; % [Hz]

% pulse repetition frequency (not provided)
prf = 250; % [Hz] not provided in paper

% sampling frequency (not used in paper)
fs = 50 * const.MHz2Hz; % [Hz]

% pulse width
pulse_width = 5e-6; % [s]

% slow time, pulse response interval (PRI)
dt_slow = 1/prf; % [s]

% create a time array corresponding to the fast time
% sampling
dt_fast_time = 1/fs;
t = (0:dt_fast_time:pulse_width-1/fs)';
range_array = t .* const.c / 2;

% define the transmitted signal
% https://en.wikipedia.org/wiki/Chirp
tx_signal = exp(2*pi*1j*((fc - B/2)*t ...
    + B/(2*pulse_width)).*t.^2);

% range and cross-range resolution
range_resolution = const.c / (2*B); % [m]
cross_range_resolution = 0.5; % [m]


%% Define the target trajectory

% define the number of scattering points off of the target
num_scatters = 1;

% define the target's initial position
target_position = [70.8, 117.5, 500];

% intialize target
target = Target(dt_slow, target_position, num_scatters);

% initialize target velocity moving along the y-axis
target.velocity = [0, 100, 0]; % [m/s]
% target.position = target.position * const.km2m;

% from the paper the target is traversing along a circular
% parth (target radius not provided)
target.radius = 30; % [m]

% calculate the period of the target's rotation
T = 1;

% calculate the number of pulses
num_cross_range_bins = round(T / dt_slow);
num_range_bins = size(t,1);


%% Propagate

% initialize received signal matrix
% [number of range cells x number of cross range cells]
rx_signal = zeros(num_range_bins, num_cross_range_bins);
target_positions = zeros(num_cross_range_bins, 2);
ranges = zeros(num_cross_range_bins, 1);

% iterate through each pulse (column)
for ipulse = 1:num_cross_range_bins

    % iterate through each target scatter
    for ipt = 1:target.num_scatters

        % calculate the range to the target center, we assuming
        % that the ground station is positioned at the origin
        range = norm(target.scatter_positions(ipt, :));
        ranges(ipulse) = range;

        % find that range index corresponding to the target
        % that minimizes the difference between each range
        % bin and the target's range
        [~, target_idx] = min(abs(range_array - range));
    
        % calculate the time delay
        tau = 2 * range / const.c;

        % solve for the phase of the signal baseband
        phase = 2 * pi * fc * tau; % need to include radial velocity
    
        % assume that only the target scattering points
        % reflect the signal and for now assume the
        % reflectivity is perfect
        reflectivity = 1;
        rx_signal(target_idx,ipulse) = ...
            reflectivity * exp(1j * phase);

    end
    
    % save target position data
    target_positions(ipulse, :) = target.position(1:2);

    % propagate target
    target.propagate();

end

%% ISAR Signal Processing
% https://www.numberanalytics.com/blog/sar-signal-processing-essentials

% range compression via match filtering using the
% transmitted signal pulse
h = conj(flipud(tx_signal));
rx_signal_range_compressed = zeros(size(rx_signal));
for irange =1:num_range_bins
    rx_signal_range_compressed(irange,:) = ...
        conv(rx_signal(irange,:), h, 'same');
end
% perform Doppler processing (or az FFT) on the range
% compressed data   
doppler_spectrum = ...
    fftshift(fft(rx_signal_range_compressed, [], 2), 2);

%% Plotting

figure
subplot(1,2,1)
plot(target_positions(:,1), target_positions(:,2))
subplot(1,2,2)
plot(ranges)

figure

subplot(1,3,1)
imagesc(abs(rx_signal));
colorbar
xlabel('cross range [m]')
ylabel('range [m]')
title('rx signal prior to compression')
axis square

subplot(1,3,2)
imagesc(abs(rx_signal_range_compressed))
colorbar
xlabel('cross range [m]')
ylabel('range [m]')
title('range compressed rx signal')
axis square

subplot(1,3,3)
imagesc(abs(doppler_spectrum))
colorbar
xlabel('Doppler [m]')
ylabel('range [m]')
title('Doppler spectrum (Range-Doppler Map)')
axis square



