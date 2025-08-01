

% instantiate plotting
p = plotting();
p.visible                       = false;
p.bool_plot_all                 = true;
p.bool_video_target_trajectory  = false;

% instantiate constants
const = Constants;

%% Define the signal model

% initialize radar parameters and LFM signal model
fc  = 10 * const.GHz2Hz; % [Hz] center frequency - X-band
B   = 149.9 * const.MHz2Hz; % [Hz] bandwidth
prf = 1 * const.kHz2Hz; % [Hz] pulse repetition frequency
fs  = 300 * const.MHz2Hz; % [Hz] sampling frequency
Tp  = 5 * const.us2s; % [s] pulse width

% in order to pad the fast time array, a rough guess at the
% maximum range must be set
max_range = 520;

% define the length of the simulation
T = 2; % [s]

% define the LFM signal
signal = LFM_Signal(fc, B, prf, fs, Tp, T, max_range);
    

%% Define the target model
    
% define the number of scattering points off of the target
num_scatters = 3;

% define the target's initial position
target_position = [70.8, 117.5, 500];

% intialize target
target = Target(signal.dt_slow, target_position, num_scatters);

% set straight line velocity towards the radar
target.velocity = [0, -10, 0];
target.p = pi/16; % [rad/s] rotate 360 degrees over the course of the simulation
target.p_dot = pi/16; % [rad/s/s]


%% Execute the imager

% execute the isar imager
output_struct = isar_imager(signal, target);

% plot
% p.plot(output_struct);









