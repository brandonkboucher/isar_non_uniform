

function [rx_signal, ranges, phases,scatterer_positions] =...
    form_raw_isar_data(...
    signal, ...
    target)

    const = Constants();
    
    %% raw isar image
    
    % initialize received signal matrix
    % [number of range cells x number of cross range cells]
    rx_signal = zeros(signal.num_pulses, signal.num_range_bins);
    
    % output vectors
    scatterer_positions = zeros(signal.num_pulses, target.num_scatterers, 3);
    ranges = zeros(signal.num_pulses, target.num_scatterers);
    phases = zeros(signal.num_pulses, target.num_scatterers);
    fprintf('Propagating simulation\n')
    
    % iterate through each pulse (column)
    for ipulse = 1:signal.num_pulses
    
        if mod(ipulse, signal.num_pulses/10) == 0
        fprintf('   %i percent complete\n', (ipulse*100/signal.num_pulses))
        end
        
        % iterate through each target scatter
        for ipt = 1:target.num_scatterers
        
            % formulate the rotation matrix about the z axis
            R = [cos(target.yaws(ipulse)), -sin(target.yaws(ipulse)), 0; ...
                 sin(target.yaws(ipulse)), cos(target.yaws(ipulse)), 0; ...
                 0, 0, 1];
            
            % calculate the scatterer's position relative to the
            % radar
            scatterer_absolute_position = ...
                (R * target.scatter_relative_positions(ipt,:)')' ...
                + target.target_center_position;
            range = norm(scatterer_absolute_position);
    
            % calculate the dealy
            delay = 2 * range / const.c; % [s]
            PRI = 1 / signal.prf;
            delay = mod(delay, PRI);

            % define the reflectivity of the scatterer
            A = 1;
            
            % define the phase 
            phase = exp( -1j * 4 * pi * signal.fc * range / const.c);
            
            % k = mod(floor(delay * signal.fs), signal.num_range_bins) + 1;
            % rx_signal(ipulse,k) = ...
            %     rx_signal(ipulse,k) + A * phase;

            % define the echo, equation 2 - not sure if this
            % should be normalized or unnormalized sinc
            scatterer_signal = ...
                A * sinc(signal.B * (signal.t_fast - delay)) * phase;
            
            % add the scatterer's contribution to the range
            % profile
            rx_signal(ipulse, :) = ...
                rx_signal(ipulse, :) + scatterer_signal.';
            
            % save target position data
            ranges(ipulse,ipt) = range;
            phases(ipulse,ipt) = phase;
    
            scatterer_position = ...
                R * target.scatter_relative_positions(ipt,:)';
            scatterer_positions(ipulse, ipt, :) = ...
                scatterer_position';
        end
    end
end

