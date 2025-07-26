function [rx_signal, output_struct] = form_raw_image( signal, target, output_struct)
    
    %% Raw ISAR data
    
    % initialize constants
    const = Constants();

    num_chirp_bins = size(signal.tx_signal, 1);

    % initialize received signal matrix
    % [number of range cells x number of cross range cells]
    rx_signal = zeros(signal.num_pulses, signal.num_range_bins);
    R0 = zeros(target.num_scatters,1);
    
    % output vectors
    scatterer_positions = zeros(signal.num_pulses, size(target.position, 2), target.num_scatters);
    los_velocities      = zeros(signal.num_pulses, target.num_scatters);
    target_positions    = zeros(signal.num_pulses, size(target.position, 2));
    ranges              = zeros(signal.num_pulses, target.num_scatters);
    
    fprintf('Propagating simulation\n')
    
    % iterate through each pulse (column)
    for ipulse = 1:signal.num_pulses
    
        if mod(ipulse, signal.num_pulses/10) == 0
            fprintf('   %i percent complete\n', (ipulse*100/signal.num_pulses))
        end
    
        % iterate through each target scatter
        for ipt = 1:target.num_scatters
    
            % calculate the range to the target center, we assuming
            % that the ground station is positioned at the origin
            % range = norm(target.scatter_positions(ipt, :));
            % ranges(ipulse) = range;
    
            % calculate the target's range
            if ipulse == 1
                R0(ipt, :) = norm(target.scatter_positions(ipt, :));
            end
            % calculate the normalized los vector
            los_vector = target.scatter_positions(ipt, :) ...
                / norm(target.scatter_positions(ipt, :));
            
            % calculate the los velocity
            los_velocities(ipulse, ipt) = ...
                dot(target.velocity, los_vector);
    
            % calculate the range resulting from a change in the
            % los velocity
            range = R0(ipt, :) ...
                + sum(los_velocities(:,ipt)) .* signal.dt_slow;
            ranges(ipulse, ipt) = range;
    
            % calculate the instantaneous Doppler frequency from (8)
            % https://en.wikipedia.org/wiki/Doppler_effect
            % https://digital-library.theiet.org/doi/10.1049/sbra504e_ch1
    
            % the Doppler frequency is defined as (2/lambda)*d
            % range / dt, which is the velocity in the radial
            % direction
            fd = 2 * signal.fc * los_velocities(ipulse, ipt) / const.c;
        
            % calculate the time delay
            tau = 2 * range / const.c;
    
            if tau < 0 || tau > signal.Tp
                warning('delay out of bounds')
            end
    
            % ensure Doppler shift begins where the target
            % returns a signal
            t_echo = (0:num_chirp_bins-1)' * signal.dt_fast_time + tau;
    
            % Doppler correction term
            tx_doppler = exp(1j * 2 * pi * fd * t_echo);
    
            % assume that only the target scattering points
            % reflect the signal and for now assume the
            % reflectivity is perfect
            reflectivity = 1;
    
            % determine index corresponding to 
            delay_idx = round(tau / signal.dt_fast_time);
    
            % initialize echo as zeros
            echo = zeros(signal.num_range_bins, 1);
            
            % ensure the echo is completely encapsulated by the
            % total number of range bins
            if delay_idx + num_chirp_bins - 1 <= signal.num_range_bins
    
                % the echo is the transmitted signal modified by
                % the target's reflectivity and the Doppler
                % shift resulting from the target's radial
                % motion
                echo(delay_idx : delay_idx + num_chirp_bins - 1) = ...
                    echo(delay_idx : delay_idx + num_chirp_bins - 1) ...    
                    + reflectivity * signal.tx_signal .* tx_doppler;
            else
                warning('Echo extends beyond fast time window.')
            end
            
            rx_signal(ipulse, :) = rx_signal(ipulse, :) + echo.';
    
            % model the received LFM signal shifted by the delay
            % [3]
            % rect = abs((t_hat - tau)/Tp) <= 1/2; 
    
            % the received signal originating from scatterer
            % rx_signal_scatterer = ...
            %     rect .*exp(pi * 1j * ...
            %     ( mu .* (t_hat - tau).^2 )); % baseband
    
            % continuous wave
            % rx_signal_scatterer = ...
            %     exp( 2 * pi * -1j * fc * (t_hat - tau)); % continuous wave
    
            % sum for total received signal
            % rx_signal(ipulse, :) = ...
            %     rx_signal(ipulse, :) + rx_signal_scatterer';
    
            % save target position data
            scatterer_positions(ipulse, ipt, :) = ...
                target.scatter_positions(ipt, :);
            target_positions(ipulse, :) = ...
                target.position;
            
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

