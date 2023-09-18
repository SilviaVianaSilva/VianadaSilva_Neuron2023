function signal_filtered = BandpassFilter(signal, Fs, Fpass)
Wn_theta = [Fpass(1)/(Fs/2) Fpass(2)/(Fs/2)]; % normalize by nyquist 
[btheta,atheta] = butter(3,Wn_theta);
signal_filtered = filtfilt(btheta,atheta,signal);