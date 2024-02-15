function [locInfo, maze] = mazeLocInfo(pathData,mazeType)
if nargin<2
    mazeType = 'fig8';
end
[pathData, maze] = addMazeLocs(pathData,mazeType);
locInfo = struct([]);
nPaths = length(pathData);
for path = 1:nPaths
    locID = pathData(path).locID;
    unqLoc = unique(locID);
    locLabels =pathData(path).locLabel;
    inds = zeros(0,2);
    labelSeq = [];
    IDSeq = [];
    for i = unqLoc'
        iLoc = locID == i;
       theseBnds = continuousRunsOfTrue(iLoc');
       nBnds = size(theseBnds,1);
       thisLabel = repmat(locLabels(theseBnds(1,1)),nBnds,1);
       thisID = repmat(i,nBnds,1);
       labelSeq = cat(1,labelSeq,thisLabel);
       IDSeq = cat(1,IDSeq,thisID);
       inds = cat(1,inds,theseBnds);       
    end
    [~,iSort] = sort(inds(:,1));
    labelSeq = labelSeq(iSort);
    IDSeq = IDSeq(iSort);
    inds = inds(iSort,:);
    nBnds = size(inds,1);
    iBool = false(nBnds,1);
    for b = 2:nBnds-1
        iBool(b) = strcmp(locLabels{b-1},locLabels{b+1});       
    end
	
    [labelSeq, idxnodes] = fixsuccnodeseq(fillArmHole(labelSeq));
	[fixedLabelSeq, idxarms] = fixsuccarmseq(fillSeqHole(labelSeq));
    a = ~idxnodes;
    b = find(~idxarms);
    idx = ~num2logical(b(a), length(fixedLabelSeq));
	if any(idx)
		fixedInds = zeros(length(idx), 2);
		fixedIDSeq = zeros(length(idx), 1);
		
		fixedInds(~idx, :) = inds;
		fixedIDSeq(~idx) = IDSeq;
		
		idx = find(idx);
		for ll = 1:length(idx)
			nextNonZeroIdx = idx(ll)+find(fixedInds(idx(ll):end, 1), 1)-1;
			fixedInds(idx(ll)-1, 2) = fixedInds(idx(ll)-1, 2) - (nextNonZeroIdx - idx(ll));
			fixedInds(idx(ll), :) = [fixedInds(idx(ll)-1, 2), fixedInds(nextNonZeroIdx, 1)];
		end
		LocRef = fieldnames(maze.locs);
		fixedIDSeq(idx) = cellfun(@(lbl) find(strcmp(LocRef, lbl)), fixedLabelSeq(idx));
	else
		fixedIDSeq = IDSeq;
		fixedInds = inds;
	end
	
    locInfo(path,1).IDSeq = fixedIDSeq;
    locInfo(path,1).labelSeq = fixedLabelSeq;
    locInfo(path,1).inds = fixedInds;
    locInfo(path,1).tInt = pathData(path).t(fixedInds);
    locInfo(path,1).duration = diff(pathData(path).t(fixedInds),1,2);
end