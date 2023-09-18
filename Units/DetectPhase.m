function [Phase, InstCycleFrequency, PerCycleFreq, signal_filtered] = DetectPhase(eeg_raw_ds,eeg_fs_ds)
MinFreq = 6;
MaxFreq = 14;
MPD = 1/MaxFreq*eeg_fs_ds;
signal_filtered = thetaphase.BandpassFilter(eeg_raw_ds, eeg_fs_ds, [5, 15]);
signal_phase = instPhase(signal_filtered);
[pks, Pk_locs] = findpeaks(signal_filtered, 'MINPEAKDISTANCE', round(MPD));
[troughs, Tr_locs] = findpeaks(signal_filtered.*-1, 'MINPEAKDISTANCE', round(MPD));

%%
peaks = Pk_locs;
troughs = Tr_locs;
Test = 0;
CycleDuration = [];
PerCycleFreq = [];

Phase = NaN(length(signal_filtered),1);
InstCycleFrequency = NaN(length(signal_filtered),1);

for i = 1:length(peaks)-1
    %if (peaks(i+1)-peaks(i))/500 < 1/MaxFreq || (peaks(i+1)-peaks(i))/500 > 1/MinFreq, continue, end  % skips peaks where next peak is outside range set by MaxFreq and MinFreq
    valley = troughs(find(troughs > peaks(i) & troughs < peaks(i+1)));
    if length(valley) ~= 1, continue, end % Makes sure there is one valley between the peaks
    % find zero crossings for descending zero and ascending zero
    [val1 ZeroCross270] = min(abs(signal_filtered(peaks(i):valley) - [(signal_filtered(peaks(i)) + signal_filtered(valley)) / 2]));
    ZeroCross270 = ZeroCross270+peaks(i)-1;
    [val1 ZeroCross90] = min(abs(signal_filtered(valley:peaks(i+1)) - [(signal_filtered(peaks(i+1)) + signal_filtered(valley)) / 2]));
    %[val2 ZeroCross90] = min(abs(signal(valley:peaks(i+1)))); %index closest to zero
    ZeroCross90 = ZeroCross90+valley-1 ;    
    
    ThetaCyclePhase = [];
    % peak to ZeroCross 270
    x = [peaks(i) ZeroCross270];
    if length(unique(x)) == 1,Test = Test+1; continue, end
    y = [180 270];
    xi = peaks(i):1:ZeroCross270; 
    yi = interp1(x,y,xi);
    
    ThetaCyclePhase(peaks(i)-peaks(i)+1:1:ZeroCross270-peaks(i)+1) = yi;
    
    % ZeroCross270 to trough
    x = [ZeroCross270 valley];
    if length(unique(x)) == 1,Test = Test+1; continue, end
    y = [270 360];
    xi = ZeroCross270:1:valley; 
    yi = interp1(x,y,xi);
    
    ThetaCyclePhase(ZeroCross270-peaks(i)+1:1:valley-peaks(i)+1) = yi;
    
    % trough to ZeroCross90 - not huge problem
    x = [valley ZeroCross90];
    if length(unique(x)) == 1,Test = Test+1; continue, end
    y = [0 90];
    xi = valley:1:ZeroCross90; 
    yi = interp1(x,y,xi);
    
    ThetaCyclePhase(valley-peaks(i)+1:1:ZeroCross90-peaks(i)+1) = yi;
    
    % ZeroCross90 to peak
    x = [ZeroCross90 peaks(i+1)];
    if length(unique(x)) == 1, Test = Test+1;continue, end
    y = [90 180];
    xi = ZeroCross90:1:peaks(i+1); 
    yi = interp1(x,y,xi);
    
    ThetaCyclePhase(ZeroCross90-peaks(i)+1:1:peaks(i+1)-peaks(i)+1) = yi;    
    
    Phase(peaks(i):peaks(i+1)) = ThetaCyclePhase;
    InstCycleFrequency(peaks(i):peaks(i+1)) = 1/((peaks(i+1)-peaks(i))/500);
    PerCycleFreq(i) = 1/((peaks(i+1)-peaks(i))/500);
end    
PerCycleFreq(PerCycleFreq == 0) = NaN;
% Test
% TESTING
figure
hold on
plot(eeg_raw_ds,'color', [0.75 0.75 0.75])
plot(signal_filtered)
plot(Phase*30,'-k')

for i = 1:length(peaks) 
plot(peaks(i),signal_filtered(peaks(i)), 'ro');  
end

for i = 1:length(Tr_locs) %  sum(Tr_locs<5000)
plot(troughs(i),signal_filtered(troughs(i)), 'ko');  
end

PercentNotIncluded  = sum(isnan(Phase))/length(signal_filtered)*100

figure, hist(Phase, 360);

% report histogram of cycle frequencies
