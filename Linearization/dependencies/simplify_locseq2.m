function [simp_seq, idx] = simplify_locseq2(locseq, lookupdirection_arms, lookupdirection_nodes)

if ~exist('lookupdirection_arms', 'var') || isempty(lookupdirection_arms)
	lookupdirection_arms = 'first';
end
if ~exist('lookupdirection_nodes', 'var') || isempty(lookupdirection_nodes)
	lookupdirection_nodes = 'first';
end

validatestring(lookupdirection_arms, {'first', 'last'});
validatestring(lookupdirection_nodes, {'first', 'last'});


% locseq = fixsuccarmseq(locseq);

idx_nodes = cleanloop(translatenode(locseq));
idx_nodes = idx_nodes(1, :);

nonrepeatNodes = locseq(idx_nodes);

armseq = cellfun(@namearm, nonrepeatNodes(1:end-1), nonrepeatNodes(2:end), 'un', 0);
armseq = [locseq{1}, armseq(:)']; % first in locseq is always OK and must be included

simp_seq = [armseq(:)'; nonrepeatNodes(:)'];
simp_seq = simp_seq(:);

idx_arms = cellfun(@(loc,i1,i2) find(strcmp(loc, locseq), 1), armseq, num2cell([1, idx_nodes(1:end-1)]), num2cell(idx_nodes));

idx = sort([idx_arms(:); idx_nodes(:)]);



function lblSeq = fillSeqHole(lblSeq)
idx = strcmpi(lblSeq, '');
idx(1) = false; idx(end) = false; % for the first and last indices we can never be sure what the previous location of the animal has been (unless we look at the videoData)

for i = find(idx(:)')
	lblSeq{i} = ['N', setdiff(intersect(lblSeq{i-1}, lblSeq{i+1}), 'A')];
end

function lblseq = translate(lblseq)

lblseq(cellfun(@isempty, lblseq)) = repmat({' '}, sum(cellfun(@isempty, lblseq)), 1);
lblseq = strrep(lblseq, ' ', 'Z');

lexicon.N1 = 'B';
lexicon.N2 = 'D';
lexicon.N3 = 'J';
lexicon.N4 = 'H';
lexicon.N5 = 'F';
lexicon.N6 = 'M';
lexicon.A16 = 'A';
lexicon.A12 = 'C';
lexicon.A23 = 'K';
lexicon.A25 = 'E';
lexicon.A34 = 'I';
lexicon.A45 = 'G';
lexicon.A56 = 'L';
lexicon.Z = '';

lblseq = cellfun(@(lbl) lexicon.(lbl), lblseq, 'un', 0);

lblseq = cat(2, lblseq{:});

function idx = lookupnode(simp_seq, locseq)
direction = 'last';
armOrNode = 'N';

idxOfInterest = find(~cellfun(@isempty, regexp(locseq, armOrNode, 'match')));

if strcmpi(direction, 'last')
	idx = idxOfInterest(regexp(translate(locseq(idxOfInterest)), ['(?<node>[' unique(translate(simp_seq)) '])(?!\k<node>)']));
else
	idx = idxOfInterest(regexp(translate(locseq(idxOfInterest)), ['(?<node>[' unique(translate(simp_seq)) '])(?!\k<node>)'])) + idxOfInterest(2) - 1;
	idx(end) = [];
	idx = [ismember(translate(locseq(idxOfInterest(1))), unique(translate(simp_seq))), idx];
end

function armName = namearm(nodeBehind, nodeForward)

nodeNumber(1)= str2double(regexp(nodeBehind, '\d{1}', 'match', 'once'));
nodeNumber(2) = str2double(regexp(nodeForward, '\d{1}', 'match', 'once'));

armName = ['A', num2str(min(nodeNumber)), num2str(max(nodeNumber))];


function p = pattern(what)
alphabet.arm = translatenode({'A12', 'A23', 'A34', 'A45', 'A56', 'A16', 'A25'});
alphabet.node = translatenode({'N1', 'N2', 'N3', 'N4', 'N5', 'N6'});


switch what
	case 'samenode'
		p = '';
		for i = 1:length(alphabet.node)
			thisNode = alphabet.node(i);
			otherNodes = setdiff(alphabet.node, thisNode);
			
			thisP = ['(' thisNode '((?=[' thisNode, alphabet.arm '])[^' otherNodes '])*' thisNode ')'];
			
			p = [thisP, '|', p];
		end
		p = ['(', p(1:end-1), ')'];
end


function l = cleanloop(string)

alphabet.arm = translatenode({'A12', 'A23', 'A34', 'A45', 'A56', 'A16', 'A25'});
alphabet.node = translatenode({'N1', 'N2', 'N3', 'N4', 'N5', 'N6'});
[s_loop, e_loop] = regexp(string, pattern('samenode'));
[startidx, endidx] = regexp(string, ['[', alphabet.node, ']']);
[~, I] = restr(startidx, s_loop+1, e_loop);
l = [startidx(~I);endidx(~I)];