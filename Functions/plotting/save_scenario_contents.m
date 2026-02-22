function save_scenario_contents(...
    sc, ...
    plotting_dir)

    if ~isfolder([plotting_dir, '/parameters'])
        mkdir([plotting_dir, '/parameters']);
    end

    %% scenario
    sc_m = ['scenarios/', class(sc), '.m'];
    save_txt = [plotting_dir, '/parameters/', class(sc), '.txt'];

    sc_content = fileread(sc_m);

    fid = fopen(save_txt, 'w');
    fprintf(fid, '%s', sc_content);
    fclose(fid);

    %% target
    target_m = ['classes/Targets/', class(sc.target), '.m'];
    target_txt = [plotting_dir, '/parameters/', class(sc.target), '.txt'];

    target_content = fileread(target_m);

    fid = fopen(target_txt, 'w');
    fprintf(fid, '%s', target_content);
    fclose(fid);

    %% run script
    run_m = 'run_isar_imager.m';
    run_txt = [plotting_dir, '/parameters/', 'run_isar_imager.txt'];

    run_content = fileread(run_m);

    fid = fopen(run_txt, 'w');
    fprintf(fid, '%s', run_content);
    fclose(fid);

end

