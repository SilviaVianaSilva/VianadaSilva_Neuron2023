function [ploc_l,ploc_r,peaks_l, peaks_r, pf_size_l, pf_size_r, splitter_cell, splitter_matrix, splitter_indx] = get_pf(cell_rate_matrix_left_in,cell_rate_matrix_right_in, min_peak_firing, min_pf_size, draw_figures)
%Determine Place Fields from Linearized Figure 8 maze, with rate Matrix from alternating
%left and right trials as input

%---------------------------------------------------------
% Matthias Haberl, 2019
%---------------------------------------------------------

% Removal of overlapping peaks between left and right trial in the stem:
% check how much percent of the peaks are overlapping, if >25% of the peak (halfwidth) in left trials is also in right trials (halfwidth), consider the same place field
splitter_cell = 0;splitter_side=0;
%% Define spliters based on bins with at least 5HZ difference 
%%maze_segment = {[1:36],[37:40],[41:60],[61:66],[67:76]};  return, delay, stem, choice, reward
use_bins = [37:61]; 

a = (cell_rate_matrix_left_in(1,use_bins)>5) .* ((cell_rate_matrix_right_in(1,use_bins)./cell_rate_matrix_left_in(1,use_bins))<0.5);  % at least 5 Hz  and at least twice as high frequ as the other turn (opposite turn <0.5 that frequ, since opposite turn often can have 0 firing)
b = (cell_rate_matrix_right_in(1,use_bins)>5) .* ((cell_rate_matrix_left_in(1,use_bins)./cell_rate_matrix_right_in(1,use_bins))<0.5);


%% splitter index calculation
% the splitter index is the inverse of a pearsons correlation between the average rate maps of left and right turns on the central arm (use_bins)
% min splitter value (= best pearsons correl of 1) is 0; max splitter value is infinite
turn_cross_cor = 1+ corrcoef((cell_rate_matrix_left_in(use_bins)), (cell_rate_matrix_right_in(use_bins)));   %adding 1 because pearson's goes from 1 to -1 (negative being inverse correlations)

%turn_cross_cor = do_correlate(cell_rate_matrix_left_in, cell_rate_matrix_right_in, use_bins);
splitter_indx = (1/turn_cross_cor(1,2)) - 0.5;  % minus 0.5 to make a pearsons corr of 1 a splitter idx of 0

rate_difference =  sum(cell_rate_matrix_left_in(1,use_bins) - cell_rate_matrix_right_in(1,use_bins));
aa = (diff([find(a)])==1); % are there neighboring bins that fulfill this criteria
bb = (diff([find(b)])==1);
splitter_cell = sum(aa)>0 | sum(bb)>0; %if on both sides no consecutive bins fulfill criterea then no splitter cell
splitter_side = sum(aa) - sum(bb); % positive value more bins on left turns fulfill criteria ("left weighted cell"), negative value more bins on right turns fullfill criteria ("right weighted cell")
splitter_weight = sum(aa) + sum(bb); % the higher this , the more unique firing for a turn
splitter_deterministics = max(sum(aa), sum(bb)); % the weight only for the side the cell is turning to/ more deterministic
splitter_matrix = [splitter_cell, splitter_side, splitter_weight, splitter_deterministics];
%% Increase rate matrix dimensions by interpolation to find the 20% values of the firing rate easier
cell_rate_matrix_left = imresize(cell_rate_matrix_left_in, [1, 100* size(cell_rate_matrix_left_in,2)],'bilinear');
cell_rate_matrix_right = imresize(cell_rate_matrix_right_in, [1, 100* size(cell_rate_matrix_right_in,2)],'bilinear');

%% All Peaks:
[peaks_l,ploc_l] = findpeaks(cell_rate_matrix_left, 'MinPeakHeight', min_peak_firing, 'SortStr', 'descend');  if isempty(peaks_l), peaks_l = 0; end
[peaks_r,ploc_r] = findpeaks(cell_rate_matrix_right, 'MinPeakHeight', min_peak_firing, 'SortStr', 'descend'); if isempty(peaks_r), peaks_r = 0; end

%% If both are empty fill with NaNs and return
if (isempty(ploc_l)&isempty(ploc_r))==1 % keep running unless both peaks are empty
[peaks_l, ploc_l, begin_pf_l, end_pf_l, pf_size_l] =   deal(NaN, NaN, NaN, NaN,NaN);
[peaks_r, ploc_r, begin_pf_r, end_pf_r, pf_size_r] =   deal(NaN, NaN, NaN, NaN,NaN);
splitter_PFs = NaN;
    return
end
%% find the dimension of each peak (20%) and calculate PF size first left then right
%%left side
if sum(peaks_l)>0 
remove_peak = zeros([numel(peaks_l),1]);
for i = 1: numel(peaks_l)
    %initialize
    end_pf_l(i) = ploc_l(i);
    begin_pf_l(i) = ploc_l(i);        
    while (cell_rate_matrix_left(end_pf_l(i))) > (0.2*  peaks_l(i)) && (end_pf_l(i) < (size(cell_rate_matrix_left,2)-10))
        end_pf_l(i) = end_pf_l(i) + 10;
    end
    while (cell_rate_matrix_left(begin_pf_l(i))) > 0.2*peaks_l(i) && begin_pf_l(i)>10
        begin_pf_l(i) = begin_pf_l(i) - 10;
    end 
end

%% Problem with this method of merging peaks, instead merge peaks if they overlap, as below
%% Merge peaks if they overlap somehow

    %{
    e_in_area = end_pf_l(i) <=(end_pf_l(1:i-1)+10) & end_pf_l(i) >= (begin_pf_l(1:i-1)-10); %10 is stepsize, need to account fo
    b_in_area = (begin_pf_l(i) <= end_pf_l(1:i-1)+10) & (begin_pf_l(i) >= begin_pf_l(1:i-1)-10); %10 is stepsize, need to account fo
    if sum(b_in_area)>0   % if any of the previous PFs covers this peak, then remove
        %disp('Merging neighboring PFs');
        for qq = find(b_in_area) %find(e_in_area) is the peak it gets merged into
            remove_peak(i) = 1; %remove this peak and add it to the larger one
            %merge PF by resetting begin to the smaller and end to the
            %larger
            begin_pf_l(qq) = min(begin_pf_l(qq), begin_pf_l(i));
            end_pf_l(qq) =  max(end_pf_l(qq), end_pf_l(i));
        end
    end
    if sum(e_in_area)>0
        for qq = find(e_in_area) %find(e_in_area) is the peak it gets merged into
            remove_peak(i) = 1;
            begin_pf_l(qq) = min(begin_pf_l(qq), begin_pf_l(i));
            end_pf_l(qq) =  max(end_pf_l(qq), end_pf_l(i));
        end
    end
%}
else
    [peaks_l, ploc_l, begin_pf_l, end_pf_l, pf_size_l] =   deal(NaN, NaN, NaN, NaN,NaN);
end

if numel(peaks_l)>1
%% remove double peaks / overlapping on left arm
[remove_peak,begin_pf_l,end_pf_l]  = remove_overlapping_pfs_within(begin_pf_l,end_pf_l);
[peaks_l, ploc_l, begin_pf_l, end_pf_l] =  deal(peaks_l(~remove_peak), ploc_l(~remove_peak), begin_pf_l(~remove_peak), end_pf_l(~remove_peak))  ;
end

pf_size_l = (2.5* (end_pf_l - begin_pf_l)) / 100; %multiply by bin size, divide by upscaling factor, to get to cm


%%right side
if sum(peaks_r)>0 
remove_peak = zeros([numel(peaks_r),1]);
for i = 1: numel(peaks_r)
    end_pf_r(i) = ploc_r(i);
    begin_pf_r(i) = ploc_r(i);
     %clear peaks that are lower (since sorted by height), if they are
    %within a already calculated PF (no need to calculate PF again)
    while cell_rate_matrix_right(end_pf_r(i)) > 0.2*peaks_r(i) && end_pf_r(i) < (size(cell_rate_matrix_right,2)-10)
        end_pf_r(i) = end_pf_r(i) + 10;
    end
    while cell_rate_matrix_right(begin_pf_r(i)) > 0.2*peaks_r(i) && begin_pf_r(i)>10
        begin_pf_r(i) = begin_pf_r(i) - 10;
    end     
    
end
else
    [peaks_r, ploc_r, begin_pf_r, end_pf_r, pf_size_r] =   deal(NaN, NaN, NaN, NaN,NaN);
end

if numel(peaks_r)>1
%% remove double peaks / overlapping on left arm
[remove_peak,begin_pf_r,end_pf_r]  = remove_overlapping_pfs_within(begin_pf_r,end_pf_r);
%remove the merged peaks
[peaks_r, ploc_r, begin_pf_r, end_pf_r] =  deal(peaks_r(~remove_peak), ploc_r(~remove_peak), begin_pf_r(~remove_peak), end_pf_r(~remove_peak))  ;
end
pf_size_r = (2.5* (end_pf_r - begin_pf_r)) / 100; 
%% Remove PF below a certain dimension
min_pf_size = min_pf_size; % in cm
keep_r = (pf_size_r>min_pf_size);
keep_l = (pf_size_l>min_pf_size);
[peaks_r, ploc_r, begin_pf_r, end_pf_r, pf_size_r] =  deal(peaks_r(keep_r), ploc_r(keep_r), begin_pf_r(keep_r), end_pf_r(keep_r),pf_size_r(keep_r));
[peaks_l, ploc_l, begin_pf_l, end_pf_l, pf_size_l] =  deal(peaks_l(keep_l), ploc_l(keep_l), begin_pf_l(keep_l), end_pf_l(keep_l),pf_size_l(keep_l));

splitter_PFs = 0; 
if (sum(peaks_r)>0)&(sum(peaks_l)>0)  %only look for PFs if both sides have a peaks
overlapping_pf = 0;
%% Check for overlap in peaks in the delay + stem
overlap_area = [37:60]; % central arm, excluded choice bifurcation; use where the mouse runs in left and right turns
l_mem = (ploc_l > min(100*overlap_area) &  ploc_l < max(100*overlap_area)); % peaks in that range
r_mem = (ploc_r > min(100*overlap_area) &  ploc_r < max(100*overlap_area)); % peaks in that range 

beg_l = round(begin_pf_l/100); beg_r = round(begin_pf_r/100); % scale back to bins for calculating overlap 
end_l = round(end_pf_l/100); end_r = round(end_pf_r/100);


if (sum(l_mem)>=1)&(sum(r_mem)>=1) % only cross ref if on both sides at least one PF in the stem area
matching_pf_rl = zeros(size(r_mem));
matching_pf_lr = zeros(size(l_mem));
keep_r = logical(ones(size(peaks_r))); keep_l = logical(ones(size(peaks_l)));
for lll = find(l_mem) %cross correlate each Pf in left trial with each Pf in right trial which had their peak in the overlap_area
    for rrr = find(r_mem)
        overl_peaks_LR = ismember([beg_l(lll):end_l(lll)], [beg_r(rrr):end_r(rrr)]);
        if sum(overl_peaks_LR)>0.25*size([beg_l(lll):end_l(lll)],2)  % check how much percent of the peaks are overlapping, if >25% of the peak (halfwidth) in left trials is also in right trials , consider the same place field
            overlapping_pf = overlapping_pf+1;
            matching_pf_rl(rrr) = lll;
            matching_pf_lr(lll) = rrr;
            if peaks_r(rrr)>=peaks_l(lll)
                keep_r(rrr) = 1;
                keep_l(lll) = 0; %discard the peak on the side that is smaller
            elseif peaks_r(rrr)<peaks_l(lll)          
                keep_r(rrr) = 0;
                keep_l(lll) = 1; %discard the peak on the side that is smaller
            end            
        end
    end
end

%now discard peaks that are calculated for both turns in the smaller side
if sum(keep_r)>0
[peaks_r, ploc_r, begin_pf_r, end_pf_r, pf_size_r] =  deal(peaks_r(keep_r), ploc_r(keep_r), begin_pf_r(keep_r), end_pf_r(keep_r),pf_size_r(keep_r));
else
[peaks_r, ploc_r, begin_pf_r, end_pf_r, pf_size_r] =  deal(NaN, NaN, NaN,NaN,NaN);   
end
if sum(keep_l)>0
[peaks_l, ploc_l, begin_pf_l, end_pf_l, pf_size_l] =  deal(peaks_l(keep_l), ploc_l(keep_l), begin_pf_l(keep_l), end_pf_l(keep_l),pf_size_l(keep_l));
else
[peaks_l, ploc_l, begin_pf_l, end_pf_l, pf_size_l] = deal(NaN, NaN, NaN,NaN,NaN); 
end

end % end if-statement requiring PF in STEM of both sides   
end % end if-statement requiring PF in both sides
%{
% splitter_PFs definition is not stringent enough, too many cells have
% non-overlapping PFs in the stem that might not be
splitter_PFs = sum(matching_pf_rl==0) + sum(matching_pf_lr==0); %check which PF left or right don't have a PF in the other turn
% splitter as defined as a cell with a PF that is only active in specific
% side turn not in the other side
if splitter_PFs == 0, splitter_side = 0,
elseif sum(matching_pf_rl==0) > sum(matching_pf_lr==0), splitter_side = 2  %right splitter cell
elseif sum(matching_pf_rl==0) <= sum(matching_pf_lr==0), splitter_side = 1 %left splitter cell
end

elseif (sum(l_mem)==0)&(sum(r_mem)>=1) 
   %PFs in right turn only
   splitter_PFs = sum(r_mem);
   splitter_side = 2;
elseif (sum(l_mem)>=1)&(sum(r_mem)==0)
   %PFs in left turn only 
   splitter_PFs = sum(l_mem);
   splitter_side = 1;

%% splitter cells
if splitter_PFs~=0, splitter_cell = 1; else, splitter_cell = 0; end 
%}
%% Total tumber of peaks
pf_num_peaks = sum(peaks_l>0) + sum(peaks_r>0);   % count number of peaks, except the 0 from the empty


%{
%% sanity check for calculation of place field size
subplot(size(avg_rate_matrix_left,1),1,ccc)
if pf_num_peaks(pf_cell_numb,1) > 0
plot(upscale_rate_matrix), hold on
else  % if no PF are found on either side plot both sides firing rate
    plot(imresize(avg_rate_matrix_left(ccc,:), [1, 100* size(avg_rate_matrix_left,2)]),'b--'),hold on
    plot(imresize(avg_rate_matrix_right(ccc,:), [1, 100* size(avg_rate_matrix_right,2)]),'k'),hold on
end
plot(pf_marker1(pf_cell_numb,1),0.2*pf_height(pf_cell_numb,1),'ro','MarkerSize',12)
plot(pf_marker2(pf_cell_numb,1),0.2*pf_height(pf_cell_numb,1),'go','MarkerSize',12)        
%}

if draw_figures == 1

num_of_figures =  length(findobj('type','figure'));
if num_of_figures < 20
figure
plot(cell_rate_matrix_left,'b--', 'LineWidth',1),hold on
plot(cell_rate_matrix_right,'r:', 'LineWidth',1),hold on

if ~isnan(ploc_r)
plot(ploc_r,peaks_r,'k*','MarkerSize',12)
plot(begin_pf_r,cell_rate_matrix_right(round(begin_pf_r)),'ro','MarkerSize',12)
plot(end_pf_r,cell_rate_matrix_right(round(end_pf_r)),'ro','MarkerSize',12)     
end
if ~isnan(ploc_l)
plot(ploc_l,peaks_l,'k*','MarkerSize',12)
plot(begin_pf_l,cell_rate_matrix_left(round(begin_pf_l)),'go','MarkerSize',12)
plot(end_pf_l,cell_rate_matrix_left(round(end_pf_l)),'go','MarkerSize',12) 
end
else
    wait(30)
end

end
%{
plot(cell_rate_matrix_left_in,'c--'),hold on
plot(cell_rate_matrix_right_in,'m:'),hold on
plot(ploc_r/100,peaks_r,'k*','MarkerSize',12)
plot(begin_pf_r/100,0.2*peaks_r,'ro','MarkerSize',12)
plot(end_pf_r/100,0.2*peaks_r,'ro','MarkerSize',12)        
plot(ploc_l/100,peaks_l,'k*','MarkerSize',12)
plot(begin_pf_l/100,0.2*peaks_l,'go','MarkerSize',12)
plot(end_pf_l/100,0.2*peaks_l,'go','MarkerSize',12)    
%}


if sum(peaks_l>0) == 0
    % Fill Nans here
        [peaks_l, ploc_l, begin_pf_l, end_pf_l, pf_size_l] =   deal(NaN, NaN, NaN, NaN,NaN);
end
if sum(peaks_r>0) == 0
    % Fill Nans here
        [peaks_r, ploc_r, begin_pf_r, end_pf_r, pf_size_r] =   deal(NaN, NaN, NaN, NaN,NaN);
end




    function [remove_peak,begin_pf,end_pf] =  remove_overlapping_pfs_within(begin_pf,end_pf )
        
        d_beg_pf = round(begin_pf/100); % scale back to bins for calculating overlap quickly
        d_end_pf = round(end_pf/100);
        remove_peak = zeros(size(begin_pf));
        for i = 1:numel(begin_pf)-1 %compare 1st to secondlast with every next until the last peak
            for pp = i+1:numel(begin_pf) %cross correlate each Pf in left trial with each Pf in right trial which had their peak in the overlap_area
                overl_peaks = ismember([d_beg_pf(i):d_end_pf(i)], [d_beg_pf(pp):d_end_pf(pp)]); %so many bins are overlapping
                if sum(overl_peaks)>0.25*size([d_beg_pf(pp):d_end_pf(pp)],2)  % check how much percent of the peaks are overlapping, if >25% of the peak (halfwidth) in left trials is also in right trials , consider the same place field
                    remove_peak(pp) = 1;
                    %replace the more accurate values now
                    begin_pf(i) = min(begin_pf(pp), begin_pf(i));
                    end_pf(i) = max(end_pf(pp), end_pf(i));
                else
                   
                end
            end
        end
    end

end
    %clear peaks that are lower (since sorted by height), if they are
    %within a already calculated PF (no need to calculate PF again) 
    % first metho for PFs removing duplicates 
    %{
    if i>1
        if sum((end_pf_l(1:i-1) > ploc_l(i)) & (ploc_l(i)> begin_pf_l(1:i-1)))>0   % if any of the previous PFs covers this peak, then remove
        remove_peak(i) = 1;
        %disp('Removing duplicate PF');        
        continue %don't bother to calculate actual PF begin and end
        end        
    end    
   
    %}
    
       function [average_corr] = do_correlate(rate_matrix_a, rate_matrix_b, areacode)
        [cells_a, pixels_a, trials_a] = size(rate_matrix_a);
        [cells_b, pixels_b, trials_b] = size(rate_matrix_b);
        %this can be done quicker  newrates = rate_matrix(:,areacode,:),
        %but errors if rate matrix is empty, would need to account for
        for cc = 1:cells_a
            for tt = 1:trials_a
                RM_a(cc,:,tt) = rate_matrix_a(cc,areacode,tt); %to take out the delay site (is 66 bins long)
            end
        end
        for cc = 1:cells_b
            for tt = 1:trials_b
                RM_b(cc,:,tt) = rate_matrix_b(cc,areacode,tt); %to take out the delay site (is 66 bins long)
            end
        end
        if cells_a ~= cells_b
            error('Mismatch in number of cells');
        end
        average_corr = [];
        for cc=1:cells_a
            corr_matrix = nan([trials_a, trials_b]); % makes a NaN matrix to be filled out
            for ta=1:trials_a
                for tb=1:trials_b
                    if ta == tb
                        corr_matrix(ta,tb) = nan; %to put the diagonal in blue
                    else
                        corr_matrix(ta,tb)= nancorr(RM_a',RM_b'); %sees if rate of change is similar in trials. if rate changes the same = 1.
                        %checkmax = nanmax([newrates(x,:,i) newrates(x,:,ii)]);
                        %if checkmax <= 0.25  matrix(i,ii) = NaN; end
                    end
                end
            end
            
            average_corr = [average_corr; nanmean(corr_matrix(:))];
        end
    end