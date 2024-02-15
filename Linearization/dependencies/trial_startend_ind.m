function start_end_ind = trial_startend_ind(tinfo)
%TRIAL_STARTEND_IND Find indices of start and end of each trial from
%trialInfo.mat
%
% start_end_ind = TRIAL_STARTEND_IND(tinfo)
%
% Note: ingores degeneracy of trials; if not desired, must be removed
% manually.

inds = tinfo.inds;
simp_ind = tinfo.simp_ind;

start_end_ind = cell2mat(cellfun(@(simp, ind) [ind(1, simp(1)), ind(2, simp(end))], simp_ind, inds, 'un', 0));