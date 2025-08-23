

% To do:
% 1. add Hamming window prior to RD to mitigate sidelobs
% 2. implement target super class

% Okay so the current issue with the ISAR imager for targets
% on the ground rotating is the selected reference range 
% profile for the range tracking algorithm is a poor 
% selection. The three scatterers are not well
% differentiated and there are a couple of artifacts. I
% think I need a more robust range tracking algorithm as
% well as motion compensation algorithm. I am also still not
% convinced that the cross range resolution of my target
% relates to the rotation extent of the target. To
% summarize, 
% 1. Range tracking: select a better reference range profile
% 2. Better motion compensation
% 3. Better test the relationship between the cross-range
% resolution of the target and the angular extent of the
% target.

% For right now, I can test the relationship between the
% cross-range resolution and the angular extent by
% simplifying the target definition to one target. I am
% fairly convinced the issue with the RD image formation is
% the phase alignment. The Doppler phases should be skewed
% with each pulse then after applying phase alignment, the
% Doppler frequency should be constant. There's a simple
% case to test whether the RD image formulation is working -
% constant LOS velocity should lead to a constant Doppler
% phase and a constant Doppler phase should give a clear
% point in the RD image. 

% I need to try rotating the targets at different
% frequencies which should equate to different Doppler
% frequencies, or at least more pronounces Doppler frequency
% differences.

% I do see cross range resolution improvement when
% increasing my simulation duration (dwell time) from 0.1
% seconds (corresponding to an angular extent of 18 deg) to
% 0.5 seconds (corresponding to an angular extent of 90
% deg). The cross-range resolution improves from 20 Hz to 10
% Hz; however, I am still skepitical of the general
% dependency on my ISAR imager's cross-range resolution and
% the angular extent of the target.

% 

% instantiate plotting
p = plotting();
p.visible                       = true;
p.bool_plot_range               = true;
p.bool_plot_range_compressed    = true;
p.bool_plot_rd                  = true;
p.bool_plot_target_trajectory2D = true;

% simulation parameters
sim_params = Simulation_Parameters();
sim_params.translational_moco   = false;
sim_params.range_doppler        = true;

% define the length of the simulation
T = 1; % [s]

% define the scenario
% sc = simple_circ_target1_sc(T);
sc = simple_circ_target1_sc(T);

% execute the isar imager
output_struct = isar_imager(sc.signal, sc.target, sim_params);

% plot the results
p.plot(output_struct);









