
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

% range and cross-range resolution
range_resolution = 1; % [m]
cross_range_resolution = 0.5; % [m]

% pulse repetition frequency (not provided)
prf = 1 * const.kHz2Hz; % [Hz] not provided in paper

% sampling frequency (not used in paper)
fs = 600 * const.MHz2Hz; % [Hz]

% slow time, pulse response interval (PRI)
dt = 1/prf; % [s]

% pulse width
pulse_width = 2e-6; % [s]

% the sampling frequency must be => 2*bandwidth according to
% the Nyquist criterion

% Configure the linear frequency modulation (LFM)
tx_signal = phased.LinearFMWaveform( ...
    'SampleRate', fs, ...
    'PulseWidth', pulse_width, ...
    'PRF', prf, ...
    'SweepBandwidth', B, ...
    'NumPulses', 1);



%% Define the target trajectory

% define the number of scattering points off of the target
num_scatters = 1;

% intialize target
target = Target(dt, num_scatters);

% initialize target velocity moving along the y-axis
target.velocity = [0, 10, 0]; % [m/s]
target.position = [70.8, 117.5, 930.2]; % [km] (from paper)
target.position = target.position * const.km2m;

% from the paper the target is traversing along a circular
% parth (target radius not provided)
target.radius = 30; % [m]

% calculate the period of the target's rotation
T = target.period;

% calculate the number of pulses
num_pulses = round(T) / dt;


%% Propagate
for ipulse = 1:num_pulses

    






end








