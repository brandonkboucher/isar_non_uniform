
%% Plotting

% instantiate plotting
p = plotting();
p.visible                       = false;
p.bool_plot_range               = true;
p.bool_plot_range_compressed    = true;
p.bool_plot_rd                  = true;
p.bool_plot_target_trajectory2D = true;

%% Simulation Parameters

% simulation parameters
sim_params = Simulation_Parameters();
sim_params.translational_moco   = false;
sim_params.range_doppler        = true;

% define the length of the simulation
T = 0.5; % [s]

%% Target

% define the scenario
sc = simple_circ_target1_sc(T);

%% Run Simulation

% execute the isar imager
output_struct = isar_imager(sc.signal, sc.target, sim_params);

% plot the results
p.plot(output_struct);









