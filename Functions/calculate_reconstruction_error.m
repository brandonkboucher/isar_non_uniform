function [err,pairs,missed,false_alarms,d] = calculate_reconstruction_error(...
    true_locations, ...
    estimated_locations, ...
    miss_penalty ... % (m) distance charged to each unmatched true scatterer
    )
%   MISS_PENALTY makes the single reported number honest when an algorithm
%   returns fewer scatterers than there are targets. Without it the RMS is
%   taken over the matched pairs alone, so an algorithm that simply declines
%   to report a hard target scores better than one that reports it badly --
%   PROMP returning one estimate for two scatterers scored 0.0034 m. Charging
%   each missed target the crossrange width of the imaged scene keeps the
%   metric a single number while pricing the miss at "could be anywhere in
%   the scene". Leave it empty for the old behaviour (missed targets ignored).

    % number of scatterers
    K = size(true_locations, 1);

    % make sure the array is in the expected dimensions
    if size(estimated_locations, 1) == 2 ...
            && size(estimated_locations, 2) ~= 2

        estimated_locations = estimated_locations.';

    elseif size(estimated_locations, 1) ~= 2 ...
            && size(estimated_locations, 2) ~= 2

        error('estimate location dimension is not (x,y)')

    end

    % euclidean distance matrix, cost matrix. computed directly rather than
    % with pdist2 so that scoring does not depend on the Statistics toolbox
    % being licensed.
    dx = true_locations(:,1) - estimated_locations(:,1).';
    dy = true_locations(:,2) - estimated_locations(:,2).';
    C = sqrt(dx.^2 + dy.^2);

    % find the matches using MATLAB's matchpairs function
    [pairs, missed, false_alarms] = matchpairs(C, 1e10);

    % calculate the total error
    d = C(sub2ind(size(C), pairs(:,1), pairs(:,2)));

    if nargin < 3 || isempty(miss_penalty)
        err = sqrt(mean(d.^2));
    else
        % every true scatterer contributes: a matched one its pairing
        % distance, an unmatched one the penalty
        n_missed = numel(missed);
        err = sqrt( (sum(d.^2) + n_missed * miss_penalty^2) ...
            / (numel(d) + n_missed) );
    end
end

