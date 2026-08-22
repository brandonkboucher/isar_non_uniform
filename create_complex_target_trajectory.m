function yaws = create_complex_target_trajectory(...
    initial_yaw,...
    initial_yawing_rate, ...
    initial_yawing_acc, ...
    yawing_jerk_magnitude, ...
    t_m)

    % define the yawing jerk
    yawing_jerk_array = yawing_jerk_magnitude * [1, -1, -1, -1];
    num_partitions = size(yawing_jerk_array,2);

    % partition into four segments
    idx_segment_incr = ...
        round(size(t_m,1)/num_partitions);
    
    yaws = zeros(size(t_m,1),1);
    yawing_rate = zeros(size(t_m,1),1);
    yawing_acc = zeros(size(t_m,1),1);
    yawing_jerk = zeros(size(t_m,1),1);

    for iseg = 1:num_partitions

        % define the segment indices
        idx = 1 + (iseg-1)*idx_segment_incr...
            :iseg*idx_segment_incr;

        if iseg == num_partitions
            if idx(end) > size(t_m,1)
                idx = idx(idx <= size(t_m,1));
            
            elseif idx(end) < size(t_m,1)
                missing_idx = 1:(size(t_m,1) - idx(end));
                idx = [idx, idx(end) + missing_idx];
            end
        end

        t_seg = t_m(idx) - t_m(idx(1));

        yaws(idx) = ...
            initial_yaw ...
            + initial_yawing_rate * t_seg ...
            + (1/2) * initial_yawing_acc * t_seg .* t_seg ...
            + (1/6) * yawing_jerk_array(iseg) * t_seg .* t_seg .* t_seg;

        yawing_rate(idx) = ...
            initial_yawing_rate ...
            + initial_yawing_acc * t_seg ...
            + (1/2) * yawing_jerk_array(iseg) * t_seg .* t_seg;

        yawing_acc(idx) = ...
            initial_yawing_acc ...
            + yawing_jerk_array(iseg) * t_seg;
        
        yawing_jerk(idx) = repmat(yawing_jerk_array(iseg), size(idx,1),1);
        
        initial_yawing_acc  = yawing_acc(idx(end));
        initial_yawing_rate = yawing_rate(idx(end));
        initial_yaw         = yaws(idx(end));
        
    end
        
end

