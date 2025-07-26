
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
clear
close all

%% Define plotting parameters

% instantiate plotting
p = plotting();
p.visible = true;
p.bool_plot_raw_isar = true;
p.bool_plot_range_compressed = true;
p.bool_plot_rd = true;
%p.bool_video_target_trajectory = true;

%% Define the tx signal parameters

% instantiate constants
const = Constants;

% center frequency - X-band
fc = 10 * const.GHz2Hz; % [Hz]

% bandwidth
B = 149.9 * const.MHz2Hz; % [Hz]

% pulse repetition frequency (not provided)
prf = 1 * const.kHz2Hz; % [Hz] this will become non-linear

% sampling frequency (not used in paper) must be greater
% than 2*Bandwidth
% https://en.wikipedia.org/wiki/Nyquist–Shannon_sampling_theorem
fs = 300 * const.MHz2Hz; % [Hz]

% pulse width
Tp = 5 * const.us2s; % [s]
max_expected_range = 520; % [m]
tau_max = 2 * max_expected_range / const.c;
padding = 1*const.us2s; % [s] additional padding

% define the chirping rate
mu = B / Tp;

% slow time, pulse response interval (PRI)
dt_slow = 1/prf; % [s]

% create a time array corresponding to the fast time
% sampling
dt_fast_time = 1/fs;
t_hat = (0:dt_fast_time:(Tp + tau_max + padding)-(1/fs))';
t_chirp = (0:dt_fast_time:Tp - (1/fs))';
range_array = t_hat .* const.c / 2;

% calculate the range resolution
range_resolution = const.c / (2 * B);

%% Define simulation parameters

% define the length of the simulation
T = 10; % [s]

% calculate the number of pulses
num_pulses = round(T / dt_slow);
num_range_bins = size(t_hat,1);
num_chirp_bins = size(t_chirp, 1);

% slow time array
t_m = (0:dt_slow:T-dt_slow)';

%% Define the target trajectory

% define the number of scattering points off of the target
num_scatters = 3;

% define the target's initial position
target_position = [70.8, 117.5, 500];

% intialize target
target = Target(dt_slow, target_position, num_scatters);

% set straight line velocity towards the radar
target.velocity = [0, -10, 0];
target.p = (2*pi) / T; % rotate 360 degrees over the course of the simulation
%target.p_dot = pi / T;

%% Define the transmitted signal

% transmitted signal
rect_tx = abs(t_chirp/Tp) <= 1/2;
tx_signal = rect_tx .* exp(1j * pi * (mu * t_chirp.^2)); % baseband

% tx_signal = exp( 2 * pi * -1j * fc * t_hat); % continuous wave

% save data
output_struct.target = target;

%% Raw ISAR data

% initialize received signal matrix
% [number of range cells x number of cross range cells]
rx_signal = zeros(num_pulses, num_range_bins);
R0 = zeros(target.num_scatters,1);

% output vectors
scatterer_positions = zeros(num_pulses, size(target.position, 2), target.num_scatters);
target_positions    = zeros(num_pulses, size(target.position, 2));
ranges              = zeros(num_pulses, target.num_scatters);

fprintf('Propagating simulation\n')

% iterate through each pulse (column)
for ipulse = 1:num_pulses

    if mod(ipulse, num_pulses/10) == 0
        fprintf('   %i percent complete\n', (ipulse*100/num_pulses))
    end

    % iterate through each target scatter
    for ipt = 1:target.num_scatters

        % calculate the range to the target center, we assuming
        % that the ground station is positioned at the origin
        % range = norm(target.scatter_positions(ipt, :));
        % ranges(ipulse) = range;

        % calculate the target's range
        if ipulse == 1
            R0(ipt, :) = norm(target.scatter_positions(ipt, :));
        end
        % calculate the normalized los vector
        los_vector = target.scatter_positions(ipt, :) ...
            / norm(target.scatter_positions(ipt, :));
        
        % calculate the los velocity
        los_velocity = dot(target.velocity, los_vector);
        
        % calculate the range resulting from a change in the
        % los velocity
        range = R0(ipt, :) + sum(los_velocity) .* dt_slow;
        ranges(ipulse, ipt) = range;

        % calculate the instantaneous Doppler frequency from (8)
        % https://en.wikipedia.org/wiki/Doppler_effect
        % https://digital-library.theiet.org/doi/10.1049/sbra504e_ch1

        % the Doppler frequency is defined as (2/lambda)*d
        % range / dt, which is the velocity in the radial
        % direction
        fd = 2 * fc * los_velocity / const.c;
    
        % calculate the time delay
        tau = 2 * range / const.c;

        if tau < 0 || tau > Tp
            warning('delay out of bounds')
        end

        % ensure Doppler shift begins where the target
        % returns a signal
        t_echo = (0:num_chirp_bins-1)' * dt_fast_time + tau;

        % Doppler correction term
        tx_doppler = exp(1j * 2 * pi * fd * t_echo);

        % assume that only the target scattering points
        % reflect the signal and for now assume the
        % reflectivity is perfect
        reflectivity = 1;

        % determine index corresponding to 
        delay_idx = round(tau / dt_fast_time);

        % double check to ensure the delay index aligns with
        % the range calculation
        [~, target_idx] = min(abs(range_array - range));

        % initialize echo as zeros
        echo = zeros(num_range_bins, 1);
        
        % ensure the echo is completely encapsulated by the
        % total number of range bins
        if delay_idx + num_chirp_bins - 1 <= num_range_bins

            % the echo is the transmitted signal modified by
            % the target's reflectivity and the Doppler
            % shift resulting from the target's radial
            % motion
            echo(delay_idx : delay_idx + num_chirp_bins - 1) = ...
                echo(delay_idx : delay_idx + num_chirp_bins - 1) ...    
                + reflectivity * tx_signal .* tx_doppler;
        else
            warning('Echo extends beyond fast time window.')
        end
        
        rx_signal(ipulse, :) = rx_signal(ipulse, :) + echo.';

        % model the received LFM signal shifted by the delay
        % [3]
        % rect = abs((t_hat - tau)/Tp) <= 1/2; 

        % the received signal originating from scatterer
        % rx_signal_scatterer = ...
        %     rect .*exp(pi * 1j * ...
        %     ( mu .* (t_hat - tau).^2 )); % baseband

        % continuous wave
        % rx_signal_scatterer = ...
        %     exp( 2 * pi * -1j * fc * (t_hat - tau)); % continuous wave

        % sum for total received signal
        % rx_signal(ipulse, :) = ...
        %     rx_signal(ipulse, :) + rx_signal_scatterer';

        % save target position data
        scatterer_positions(ipulse, ipt, :) = ...
            target.scatter_positions(ipt, :);
        target_positions(ipulse, :) = ...
            target.position;
        
    end
    
    % propagate target
    target.propagate();

end

% save data
output_struct.ranges = ranges;
output_struct.t_m = t_m;
output_struct.rx_signal = rx_signal;
output_struct.range_array = range_array;
output_struct.scatterer_positions = scatterer_positions;
output_struct.target_positions = target_positions;

fprintf('Performing post-processing.\n')
%% Range compression
% https://www.numberanalytics.com/blog/sar-signal-processing-essentials

% range compression via match filtering using the
% transmitted signal pulse
h = conj(flipud(tx_signal));
rx_signal_range_compressed = zeros(size(rx_signal));
for ipulse =1:num_pulses    
    rx_signal_range_compressed(ipulse,:) = ...
        conv(rx_signal(ipulse,:), h, "same");
end

% adjust the range axis by the match filter's group delay
range_array = range_array - (Tp * const.c / 4);
output_struct.range_array = range_array;

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


%% Phase Adjustment

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

fprintf('Plotting results.\n')

% plot
p.plot(output_struct);


fprintf('Simulation complete.\n')
