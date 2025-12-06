
% The Map-Drift algorithm estimates rotational motion
% parameters by splitting up the aperture into sub-apertures
% and evaluating the cross-correlation between sub-apertures
% to determine the drift of the scatterers.

function drift_estimation = map_drift(...
    rx_signal,...
    num_sub_apertures, ...
    ranges, ...
    signal)

    debug = false;

    % calculate the range bin spacing using the sampling
    % frequency
    const = Constants();
    range_bin_spacing = const.c / (2 * signal.fs);

    % begin by splitting the raw ISAR image into
    % sub-apertures
    [num_pulses, num_ranges] = size(rx_signal);

    % check that the number of pulses is divisible by the
    % number of sub-apertures
    num_pulses_per_subaperture = round(num_pulses/num_sub_apertures);
    if mod(num_pulses, num_sub_apertures) ~= 0
        warning('Number of pulses is not divible by the number of sub-apertures.')
    end

    % define the subapertures [nPulses x nRanges x nSubApertures]
    % subaperature_images = zeros(...
    %     num_pulses_per_subaperture, ...
    %     num_ranges, ...
    %     num_sub_apertures);

    ref_image = rx_signal(1:num_pulses_per_subaperture,:);

    [~, max_idx] = max(abs(ref_image(:)));
    [ref_x, ref_y] = ind2sub(size(ref_image), max_idx);

    % define an array maintaining the change in the range
    % bin of the target btw subapertures, this define the
    % motion of the target
    range_shifts = ...
        zeros(num_sub_apertures-1,1);

    % save the indices used to evaluate the
    % cross-correlations
    xcorr_idx = zeros(num_sub_apertures,1);
    xcorr_idx(1) = round(size(ref_image,1)/2);

    for isub = 1:num_sub_apertures-1

        % determine sub-aperture index
        idx = isub*num_pulses_per_subaperture+1 ...
            :(isub+1)*num_pulses_per_subaperture;
        subaperture_image = rx_signal(idx, :);

        % for each sub-aperture find the point which 
        % maximizes the cross-correlation
        [az_shift, range_shift] = fft_cross_correlation(...
            ref_image, subaperture_image);

        xcorr_idx(isub+1) = ...
            idx(round(size(subaperture_image,1)/2));

        ref_pulse = ref_image(xcorr_idx(1),:);
        subaperture_pulse = subaperture_image(round(size(subaperture_image,1)/2),:);

        % evaluate the cross-correlation of the center
        % pulses within each subaperture
        [c, lag] = xcorr(...
            subaperture_pulse,...
            ref_pulse);

        if debug

            % find the peak points between the reference
            % image and the sub-aperture
            [~, ref_max_idx] = max(abs(ref_pulse));
            [~, sub_max_idx] = max(abs(subaperture_pulse));

            figure
            subplot(1,2,1)
            stem(lag, abs(c))
            xlabel('lag (range bins)')
            xlim([-20, 20])
            ylabel('cross-correlation')
            title(['Subaperture ', num2str(isub+1), ' cross correlation'])

            subplot(1,2,2)
            imagesc(20*log10(abs(subaperture_image)))
            title('sub-aperture')
            xlim([4000, 4030])
            colorbar


        end

        % find the lag that corresponds to a maximization of
        % the 1D cross correlation btw range profiles this
        % corresponds to the motion of the target
        [~,maxIdx]=max(c);
        range_shifts(isub) = lag(maxIdx);

    end

    % convert range bins to meters
    range_shifts = range_shifts .* range_bin_spacing;

    %% Debugging
    % calculate the range bin changes of the targets
    % relative to the reference image (index halfway in the
    % reference image, first index of xcorr_idx)
    % range_reference = ranges(xcorr_idx(1),:);
    % range_difference = ranges(xcorr_idx(2:end),:) - range_reference;
    % 
    % % convert to range bins
    % range_bin_difference = range_difference ./ signal.range_resolution;
    % 
    % % we now have the drift resulting from the rotational
    % % motion of the target - estimate the drift over the
    % % number of pulses
    % subaperture_centers = ...
    %     round(linspace(0,num_pulses, num_sub_apertures) ...
    %     + num_pulses / (2 * num_sub_apertures));
    % 
    % % perform a first order fit
    % coeffs = polyfit(subaperture_centers(1:end-1),...
    %     range_bin_difference, 1);
    
    %%

    % we now have the drift resulting from the rotational
    % motion of the target - estimate the drift over the
    % number of pulses
    subaperture_centers = ...
        round(linspace(0,num_pulses, num_sub_apertures) ...
        + num_pulses / (2 * num_sub_apertures));

    % perform a first order fit
    coeffs = polyfit(subaperture_centers(1:end-1),...
        range_shifts, 1);

    % evaluate the drift at every pulse
    pulse_idx = 0:num_pulses - 1;
    drift_estimation = polyval(coeffs,pulse_idx);
    drift_estimation = drift_estimation';

    % normalize to find the drift estimation per pulse
    % drift_estimation = drift_estimation / num_pulses;

    if debug
        f = figure('Visible','off');
        subplot(1,2,1)
        plot(ranges(1) + drift_estimation, 'LineWidth', 2, 'DisplayName', 'Estimate')
        hold on
        plot(ranges, 'LineWidth', 2, 'DisplayName', 'Actual')
        title(['Linear drift estimation for ', num2str(num_sub_apertures), ' sub-apertures'], 'FontSize', 24)
        xlabel('Pulses', 'FontSize', 16)
        ylabel('Range [m]', 'FontSize', 16)
        legend
        subplot(1,2,2)
        plot(abs(ranges - (ranges(1) + drift_estimation)), 'LineWidth', 2)
        title(['Linear drift estimation error for ', num2str(num_sub_apertures), ' sub-apertures'], 'FontSize', 24)
        xlabel('Pulses', 'FontSize', 16)
        ylabel('Range Error [m]', 'FontSize', 16)
        set(gcf, 'Position', get(0, 'Screensize'));
        saveas(f,'plots/map_drift.png')
    end
end







