

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

% define the length of the simulation
T = 2; % [s]

% in order to pad the fast time array, a rough guess at the
% maximum range must be set
max_range = 520;

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

% define the roll acceleration array
p_dot_array = (0:pi/45:2*pi)';


%% Iterate through the acceleration array

for iacc = 1:size(p_dot_array,1)

    % set the target angular acceleration
    target.p_dot = p_dot_array(iacc);

    % execute the isar imager
    output_struct = isar_imager(signal, target);

    % plot
    p.plot(output_struct, ['analysis1/acc_', num2str(rad2deg(p_dot_array(iacc)))]);

end








