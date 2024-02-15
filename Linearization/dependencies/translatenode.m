function lblseq = translatenode(lblseq)

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
