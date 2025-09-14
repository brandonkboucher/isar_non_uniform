function [rx_signal, output_struct] = form_raw_image_basic( signal, target, output_struct)
    
    %% Raw ISAR data
    
    % initialize constants
    const = Constants();
    c = const.c;

    % initialize received signal matrix
    % [number of range cells x number of cross range cells]
    rx_signal = zeros(signal.num_pulses, signal.num_range_bins);
    
    % output vectors
    scatterer_positions = zeros(signal.num_pulses, target.num_scatterers, size(target.position, 2));
    target_positions    = zeros(signal.num_pulses, size(target.position, 2));
    ranges              = zeros(signal.num_pulses, target.num_scatterers);
    fprintf('Propagating simulation\n')
    
    % iterate through each pulse (column)
    for ipulse = 1:signal.num_pulses
    
        if mod(ipulse, signal.num_pulses/10) == 0
            fprintf('   %i percent complete\n', (ipulse*100/signal.num_pulses))
        end
    
        % array corresponding the fast time
        t_fast = signal.t_fast;

        % iterate through each target scatter
        for ipt = 1:target.num_scatterers
    
            % this simulator is based off of Novel
            % Non-Uniform Rotational Motion Estimation and
            % Compensation Method and is a tabletop
            % simulation

            % extract the position of each scatterer
            % relative to the axis of rotation
            x = target.scatter_relative_positions(ipt,1); % [m]
            y = target.scatter_relative_positions(ipt,2); % [m]
            R0 = norm(target.position); % [m]
            theta = target.yaw;

            % extract the signal information
            B = signal.B; % [Hz] Bandwidth
            fc = signal.fc; % [Hz] center frequency

            % calculate the range of the scatterer relative
            % to the radar at the origin (0,0)
            range = R0 + x*sin(theta) + y*cos(theta); % [m]
            delay = 2 * range / c; % [s]

            % define the reflectivity of the scatterer
            A = 1;

            % define the phase 
            phase = exp( -1j * 4 * pi * fc * range / c);

            % define the echo, equation 2 - not sure if this
            % should be normalized or unnormalized
            scatterer_signal = ...
                A * sinc(B * (t_fast - delay)) * phase;

            % add the scatterer's contribution to the range
            % profile
            rx_signal(ipulse, :) = ...
                rx_signal(ipulse, :) + scatterer_signal';

            

            % save target position data
            scatterer_positions(ipulse, ipt, :) = ...
                target.scatter_positions(ipt, :);
           
            target_positions(ipulse, :) = target.position;
            ranges(ipulse,ipt) = range;
            
        end
        
        % propagate target
        target.propagate();
    
    end

    % save data
    output_struct.ranges = ranges;
    output_struct.rx_signal = rx_signal;
    output_struct.scatterer_positions = scatterer_positions;
    output_struct.target_positions = target_positions;

end

