
close all

%% Plotting

% instantiate plotting
p = plotting();
p.visible                       = false;
p.bool_plot_all                 = false;

%% Simulation Parameters

% simulation parameters
sim_params = Simulation_Parameters();
sim_params.translational_moco   = false;
sim_params.backprojection       = true;
sim_params.range_doppler        = false;

% define the length of the simulation
T = 1; % [s]

%% Target

% define the scenario
sc = tabletop_L_simple_sig_sc(T);
%sc = simple_circ_target1_sc(T);

%% Run Simulation

% execute the isar imager
output_struct = isar_imager(sc.signal, sc.target, sim_params);

% plot the results
p.plot(output_struct);









