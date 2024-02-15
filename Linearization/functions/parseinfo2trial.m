function parseInfo = parseinfo2trial(parseInfo, trialInfo)

start_end_ind = arrayfun(@trialInfo2simple.trial_startend_ind, trialInfo, 'un', 0);

parseInfo = cellfun(@accFunc_f, num2cell(parseInfo), start_end_ind, 'un', 0);

function parsed = accFunc_f(parseInfo, idx)

parsed = cellfun(@(idx) extract(parseInfo, idx), row2cell(idx), 'un', 0);
parsed = cat(1, parsed{:});

function x = extract(parseInfo, idx)

x = structfun(@(x) find(x(idx(1):idx(2)))+idx(1)-1, parseInfo, 'un', 0);