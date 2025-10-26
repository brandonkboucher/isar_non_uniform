
close all

show_plots = true;
save_data = false;
save_plots = false;

% define the length of the simulation
T = 1; % [s]

% define the scenario
sc = scenario_basic(T);

% execute the isar imager
output = isar_imager(sc);

%% plotting
plot_backprojection(...
    output.rx_signal_bp, ...
    output.x_array, ...
    output.y_array, ...
    show_plots, ...
    save_plots);

plot_target_range_and_traj(...
    sc.signal.t_slow, ...
    output.ranges, ...
    output.scatterer_positions, ...
    show_plots, ...
    save_plots);

%% save data
if save_data
    save('data/rx_signal_bp_sd.mat', "output.rx_signal_bp");
end






