
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

%% Define plotting parameters

% instantiate plotting
p = plotting();
% p.plot_all = true;
p.plot_all = true;

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
num_pulses = round(T / dt_slow);
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

%% Raw ISAR data

% initialize received signal matrix
% [number of range cells x number of cross range cells]
rx_signal = zeros(num_pulses, num_range_bins);
R0 = zeros(target.num_scatters,1);

% output vectors
target_positions = zeros(num_pulses, size(target.position, 2));
los_velocities = zeros(num_pulses, 1);
ranges = zeros(num_pulses, 1);
fds = zeros(num_pulses, 1);

% iterate through each pulse (column)
for ipulse = 1:num_pulses

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

% save data
output_struct.ranges = ranges;
output_struct.rx_signal = rx_signal;
output_struct.range_array = range_array;
output_struct.target_positions = target_positions;

%% Range compression
% https://www.numberanalytics.com/blog/sar-signal-processing-essentials

% range compression via match filtering using the
% transmitted signal pulse
h = conj(flipud(tx_signal));
rx_signal_range_compressed = zeros(size(rx_signal));
for ipulse =1:num_pulses
    rx_signal_range_compressed(ipulse,:) = ...
        conv(rx_signal(ipulse,:), h, 'same');
end

% save data 
output_struct.rx_signal_range_compressed = ...
    rx_signal_range_compressed;

%% Range tracking 

% Range tracking or range alignment functions to compensate
% for the translational motion, and aligning the range
% profiles. A reference range profile is selected along a
% pulse, and the cross correlation is taken for each pulse.
% The resulting delay between the reference and each pulse
% is the amount the range profile needs to be shifted in
% order to be aligned

% select a reference pulse, the center cross-range bin
ref_pulse = round(num_pulses / 2);

% extract the corresponding range profile
ref_profile = rx_signal_range_compressed(ref_pulse, :);

% iterate over each pulse and align the range profile with
% the reference range profile
rx_signal_aligned = zeros(size(rx_signal_range_compressed));
for ipulse = 1:num_pulses

    % extract the range profile
    profile = rx_signal_range_compressed(ipulse, :);

    % evaluate the cross correlation
    [corr,lag] = xcorr(profile, ref_profile);

    % find the lag that maximizes the absolute value of the
    % correlation
    [~, max_idx] = max(abs(corr));

    % shift the range profile to maximize correlation with
    % the reference profile
    rx_signal_aligned(ipulse, :) = ...
        circshift(profile, -lag(max_idx));

end

% save data
output_struct.rx_signal_aligned = rx_signal_aligned;

%% Phase Adjustment (autofocus)

% Phase adjustment is a way to sharpen the image along the
% cross-range direction prior to Range-Doppler processing,
% and works by assuming the phase of the dominant scatters
% is roughly constant. Any phase change across pulses is
% assumed error or phase drift. To correct, find the
% dominant scatters phase angles and average the phases for
% each pulse.

% I'm having a real tough time finding a coherent, fully 
% written out version of the MVM so I'm going to rely on the
% AI to help for now. If I can't get this to work then I'll
% try implementing the phase gradient algorithm

% select number of bright scatters
K = 1; % we only should have one bright scatterer

% find the index of the most dominant scatterer along the
% range axis
[~, dominant_idx] = maxk(max(abs(rx_signal_aligned), [], 1), K);

% extract the phases of the most dominant scatter across
% pulses
phases = angle(rx_signal_aligned(:,dominant_idx));

% apply the phase correction over the range bins
correction = exp(-1j * phases);
rx_autofocused = rx_signal_aligned .* correction;

% save data
output_struct.rx_autofocused = rx_autofocused;

%% RD processing

% perform Doppler processing (or az FFT) on the range
% compressed data   
rx_signal_rd = ...
    fftshift(fft(rx_autofocused, [], 1), 1);

% save data
output_struct.rx_signal_rd = rx_signal_rd;

% plot
p.plot(output_struct);
