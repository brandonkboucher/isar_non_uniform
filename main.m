
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

% [1] Implementation of ISAR Imaging with Step Frequency and LFM Waveforms Using Gabor Transform
% [2] https://en.wikipedia.org/wiki/Chirp
% [3] ISAR Imaging of Targets With Complex Motion Based on Discrete Chirp Fourier Transform for Cubic Chirps
% [4] https://ieeexplore.ieee.org/abstract/document/7985358
close all

%% Define the tx signal parameters

% instantiate plotting
p = plotting();

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
Tp = 5e-6; % [s]

% define the chirping rate
mu = B / Tp;

% slow time, pulse response interval (PRI)
dt_slow = 1/prf; % [s]

% create a time array corresponding to the fast time
% sampling
dt_fast_time = 1/fs;
t_hat = (0:dt_fast_time:Tp-1/fs)';
range_array = t_hat .* const.c / 2;

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

% set straight line velocity towards the radar
target.velocity = [-10, -10, 0];
target.acceleration = [-2,-1,0];
target.velocity_magnitude = norm(target.velocity);

% from the paper the target is traversing along a circular
% parth (target radius not provided)
% target.radius = 30; % [m]

% initialize target velocity magnitude
% target.velocity_magnitude = 100; % [m/s]

% calculate the period of the target's rotation
T = 10; % [s]

% calculate the number of pulses
num_cross_range_bins = round(T / dt_slow);
num_range_bins = size(t_hat,1);

% slow time array
t_m = (0:dt_slow:T-dt_slow)';

% full time array t = fast time + slow time [2]
t = t_hat + t_m';

% define the transmitted signal [4]
rect = abs(t_hat/Tp) <= 1/2; 
tx_signal = rect .* exp(pi * 1j * ...
    ( 2 * fc * t_hat + mu .* t_hat.^2 ));
% tx_signal = exp( 2 * pi * -1j * fc * t_hat); % continuous wave

%% Propagate

% initialize received signal matrix
% [number of range cells x number of cross range cells]
rx_signal = zeros(num_cross_range_bins, num_range_bins);
R0 = zeros(target.num_scatters,1);

% output vectors
target_positions = zeros(num_cross_range_bins, size(target.position, 2));
los_velocities = zeros(num_cross_range_bins, 1);
ranges = zeros(num_cross_range_bins, 1);
fds = zeros(num_cross_range_bins, 1);

% iterate through each pulse (column)
for ipulse = 1:num_cross_range_bins

    % iterate through each target scatter
    for ipt = 1:target.num_scatters

        % calculate the range to the target center, we assuming
        % that the ground station is positioned at the origin
        % range = norm(target.scatter_positions(ipt, :));
        % ranges(ipulse) = range;

        % explicit Doppler representation, determining the
        % range using the velocity along the line-of-sight
        if ipulse == 1
            R0(ipt, :) = norm(target.scatter_positions(ipt, :));
        end
        los_vector = target.scatter_positions(ipt, :) ...
            / norm(target.scatter_positions(ipt, :));
        los_velocities(ipulse) = ...
            dot(target.velocity, los_vector);
        range = R0(ipt, :) ...
            + sum(los_velocities) .* dt_slow;
        ranges(ipulse) = range;

        % calculate the instantaneous Doppler frequency from (8)
        % https://en.wikipedia.org/wiki/Doppler_effect
        % https://digital-library.theiet.org/doi/10.1049/sbra504e_ch1

        % the Doppler frequency is defined as (2/lambda)*d
        % range / dt, which is the velocity in the radial
        % direction
        fd = 2 * fc * los_velocities(ipulse) / const.c;
        fds(ipulse) = fd;
    
        % calculate the time delay
        tau = 2 * range / const.c;

        % assume that only the target scattering points
        % reflect the signal and for now assume the
        % reflectivity is perfect
        reflectivity = 1;

        % model the received LFM signal shifted by the delay
        % [3]
        rect = abs((t_hat - tau)/Tp) <= 1/2; 

        rx_signal(ipulse, :) = ...
            rect .* exp(pi * 1j * ...
            ( 2 * fc * (t_hat - tau) ...
            + mu .* (t_hat - tau).^2 ));

    end
    
    % save target position data
    target_positions(ipulse, :) = target.position;

    % propagate target
    target.propagate();

end

%% Propagate plots

figure
subplot(1,2,1)
imagesc(abs(rx_signal))
title('Absolute value of received signal - raw ISAR data')
xlabel('range bin')
ylabel('pulse index')
colorbar

subplot(1,2,2)
imagesc(real(rx_signal))
title('Real value of received signal - raw ISAR data')
xlabel('range bin')
ylabel('pulse index')
colorbar

figure
subplot(1,2,1)
plot(range_array, imag(rx_signal(1500, :)))
hold on
xline(ranges(1500), 'Color', 'r', 'LineWidth', 2, 'LineStyle','--')
title('Imaginary value of received signal for pulse 1500 - raw ISAR data')
xlabel('range [m]')
ylabel('Amplitude')

subplot(1,2,2)
plot(range_array, real(rx_signal(1500, :)))
hold on
xline(ranges(1500), 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--')
title('Real value of received signal for pulse 1500 - raw ISAR data')
xlabel('range [m]')
ylabel('Amplitude')

figure
s = surf(real(rx_signal));
s.EdgeColor = 'none';
xlabel('range bins')
ylabel('pulses')

p.plot_range(ranges)
p.plot_trajectory(target_positions);

%% ISAR Signal Processing
% https://www.numberanalytics.com/blog/sar-signal-processing-essentials

% range compression via match filtering using the
% transmitted signal pulse
h = conj(flipud(tx_signal));
rx_signal_range_compressed = zeros(size(rx_signal));
for ipulse =1:num_cross_range_bins
    rx_signal_range_compressed(ipulse,:) = ...
        conv(rx_signal(ipulse,:), h, 'same');
end

%% ISAR Signal Processing plots

figure
subplot(1,2,1)
imagesc(abs(rx_signal_range_compressed))
title('Absolute value range compressed ISAR image')
xlabel('range index')
ylabel('pulse index')
colorbar

subplot(1,2,2)
imagesc(real(rx_signal_range_compressed))
title('Real Range compressed ISAR image')
xlabel('range index')
ylabel('pulse index')
colorbar

%% Motion compensation

% compensate for the motion along the line-of-sight
moco_doppler_spectrum = zeros(size(doppler_spectrum));
phase_corrections = zeros(num_cross_range_bins, 1);
for ipulse = 1:num_cross_range_bins

    % calculate the instantaneous Doppler frequency from (8)
    % https://en.wikipedia.org/wiki/Doppler_effect
    % fd = 2 * fc * los_velocities(ipulse) / const.c;

    % calculate the phase resulting from the line-of-sight
    % velocity of the target
    % phase_correction = ...
    %     exp(-1j * 2 * pi * fd * (ipulse - 1) * dt_slow);
    translational_phase_correction = ...
        exp(-1j * 4 * pi * fc * ...
        (R0 + los_velocities(ipulse))/const.c);

    phase_corrections(ipulse) = ...
        angle(translational_phase_correction);

    % calculate the motion compensated spectrum
    moco_doppler_spectrum(:, ipulse) = ...
        doppler_spectrum(:,ipulse) .* translational_phase_correction;

end

%% RD processing

% perform Doppler processing (or az FFT) on the range
% compressed data   
doppler_spectrum = ...
    fftshift(fft(rx_signal_range_compressed, [], 2), 2);

fft_178 = fftshift(fft(abs(rx_signal_range_compressed(178, :))));

figure
plot(1:250, fft_178)

%% Plotting

figure
subplot(1,2,1)
plot(target_positions(:,1), target_positions(:,2))
axis square
subplot(1,2,2)
plot(ranges)
axis square

figure

subplot(2,2,1)
imagesc(real(rx_signal));
colorbar
xlabel('cross range [m]')
ylabel('range [m]')
title('rx signal prior to compression')
axis square

subplot(2,2,2)
imagesc(real(rx_signal_range_compressed))
colorbar
xlabel('cross range [m]')
ylabel('range [m]')
title('range compressed rx signal')
axis square

subplot(2,2,3)
imagesc(real(doppler_spectrum))
colorbar
xlabel('Doppler [m]')
ylabel('range [m]')
title('Doppler spectrum (Range-Doppler Map)')
axis square

subplot(2,2,4)
imagesc(real(moco_doppler_spectrum))
colorbar
xlabel('Doppler [m]')
ylabel('range [m]')
title('Motion compensated Doppler spectrum (Range-Doppler Map)')
axis square

figure
plot(abs(doppler_spectrum(:)), 'DisplayName', 'before')
hold on
plot(abs(moco_doppler_spectrum(:)), 'DisplayName', 'after')
legend

figure
plot(abs(doppler_spectrum(:)) - abs(moco_doppler_spectrum(:)), 'DisplayName', 'before')

figure
plot(1:num_cross_range_bins, phase_corrections)
xlabel('Pulse index')
ylabel('Phase (radians)')
title('Phase correction over time')

figure
subplot(1,3,1)
plot(1:num_cross_range_bins, fds)
title('doppler')
subplot(1,3,2)
plot(1:num_cross_range_bins, ranges)
title('ranges')
subplot(1,3,3)
plot(1:num_cross_range_bins, los_velocities)
title('los velocity')

target_range_idx = 178; % find the index close to target range
figure
plot(angle(rx_signal(target_range_idx, :)));
title('Phase over pulses at target range');
xlabel('Pulse index');
ylabel('Phase (radians)');