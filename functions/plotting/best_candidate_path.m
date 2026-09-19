function [pp, n_acc] = best_candidate_path(c, amb)
% BEST_CANDIDATE_PATH  Seed plus Newton iterates of the highest-correlation
% candidate in ambiguity AMB. N_ACC counts the seed and the accepted iterates;
% a final column beyond N_ACC is the step Newton rejected (the history records
% a step before the acceptance test).

    idx     = find(c.amb == amb);
    [~, j]  = max(c.corr(idx));
    ic      = idx(j);
    pp      = [c.seed(:,ic), c.path{ic}];
    n_acc   = size(pp,2) - (size(pp,2) > 1 && any(pp(:,end) ~= c.p(:,ic)));
end
