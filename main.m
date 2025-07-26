
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
p.visible = false;
p.bool_plot_all = true;
p.bool_video_target_trajectory = true;

%% Define the tx signal model

% instantiate constants
const = Constants;

% initialize radar parameters and LFM signal model
fc  = 10 * const.GHz2Hz; % [Hz] center frequency - X-band
B   = 149.9 * const.MHz2Hz; % [Hz] bandwidth
prf = 1 * const.kHz2Hz; % [Hz] pulse repetition frequency
fs  = 300 * const.MHz2Hz; % [Hz] sampling frequency
Tp  = 5 * const.us2s; % [s] pulse width

% define the length of the simulation
T = 10; % [s]

% define the LFM signal
signal = LFM_Signal(fc, B, prf, fs, Tp, T);

% slow time array
t_m = (0:signal.dt_slow:T-signal.dt_slow)';

%% Define the target model

% define the number of scattering points off of the target
num_scatters = 3;

% define the target's initial position
target_position = [70.8, 117.5, 500];

% intialize target
target = Target(signal.dt_slow, target_position, num_scatters);

% set straight line velocity towards the radar
target.velocity = [0, -10, 0];
target.p = (2*pi) / T; % rotate 360 degrees over the course of the simulation
%target.p_dot = pi / T;

% save data
output_struct.signal = signal;
output_struct.target = target;
output_struct.t_m = t_m;

%% ISAR Imaging

% formulate the raw isar image
[rx_signal, output_struct] = ...
    form_raw_image( signal, target, output_struct);

fprintf('Performing post-processing.\n')

% perform range compression via a matched filter by
% evaluating a convolution of the received signal with a
% time reverse conjugate of the transmitted signal
fprintf('   Range compression.\n')
[rx_signal_range_compressed,output_struct] = ...
    range_compression(...
    signal, ...
    rx_signal, ...
    output_struct);

% align the range profiles across all pulses relative to a
% reference range profile
fprintf('   Range alignment.\n')
[rx_signal_aligned,output_struct] = range_tracking( ...
    signal, ...
    rx_signal_range_compressed, ...
    output_struct);

% perform coarse phase correction based on dominant
% scatterers
fprintf('   Phase adjustment.\n')
[rx_autofocused,output_struct] = phase_adjustment(...
    rx_signal_aligned, ...
    1, ... % number of bright scatterers
    output_struct);

% perform Range-Doppler processing across pulses to
% formulate the final ISAR image
fprintf('   Range-Doppler.\n')
[rx_signal_rd,output_struct] = rd_processing( ...
    rx_autofocused,...
    output_struct);

%% Plot results

fprintf('Plotting results.\n')

% plot
p.plot(output_struct);


fprintf('Simulation complete.\n')
