function [freqPeak,mod_score,lag,cor] = AC_Theta(tSp, img_outfile, plot_figures)
addpath(genpath('C:\Users\Silvia\Dropbox\SilviaProjectCode\FFTtest'));

%%Initialize variables
Burst = 3;
Frequency = [9:0.005:14];
Cyclenum = 1500;
timediff = 0.020;
startSp = 1500; 

r2 = randi(Burst,Cyclenum,1);
r2(1) = 1; %% Assume the first spike in the spike train is not a burst. 
randBurst = randi(2,Cyclenum,1);

t_bin = 0.005; %% Time resolution for the autocorrelation

     [cor, lag] = CrossCorrel(tSp', tSp');     
     [cor, lag, smooth] = SmoothCor(cor, lag, t_bin); % solve for smoothed signal, and chop off center peak and negative half of autocorrelation

%end

%% Now calculate intrinsic frequency from the FFT (without the first burst)
     [S, f] = getChronuxSpectrum(cor, t_bin);
     Sall = S(f>6.2&f<14);
     fall = f(f>6.2&f<14);
     [~,freqP] = max(Sall);
     %freqPeak(currf) = fall(freqP);
     freqPeak = fall(freqP);  %intrinsic frequency
clear S f        

%% remove first 50ms (to remove the burst) then calculate
lag_rm_burst = lag(lag>0.05&lag<=0.6);
cor_rm_burst = cor(lag>0.05&lag<=0.6);

%% Now calculate  modulation score from the FFT (including the first burst)
epoch_ind =1;
     [S, f] = getChronuxSpectrum(cor_rm_burst, t_bin);  
    [xmax, imax] = extrema(S); % find local maxima and see if any lie within theta

    ftmp = f(imax); % all frequ peaks
    stmp = S(imax); % strength, all frequ peaks

    [a_peak,ind] = max(stmp(ftmp > 7 & ftmp < 14)); % highest peak within theta
    ftmp = ftmp(ftmp > 7 & ftmp < 14);

    if ~isempty(ind)                                    % if theta peak was found
        f_peak(epoch_ind) = ftmp(ind)';     
        a_peak_av = mean(S(f>f_peak(epoch_ind)-1 & f<f_peak(epoch_ind)+1));
        peak = S(f==f_peak(epoch_ind));
    else                                                % if theta peak was not found
        f_peak(epoch_ind) = 0;
        peak = 0;
        a_peak_av = 0;
    end 
    
    power_ratio(epoch_ind) = a_peak_av / mean(S);  %divide by the average FFT
    mod_score = power_ratio;
    
    if ~(mod_score>=0)  % if no spikes after first peak in autocorrelation the spectogram is 0 and mod_score will become NaN -> set to 0 since cell is not modulated
    mod_score = 0;   
    end
    
%%    Plot and add info
     if plot_figures == 1
     fig1 = figure()     
     plot(lag,cor, 'k', 'LineWidth',2)
     %print(fig1, '-dpng', strcat(img_outfile,'.png'));
     img_outfile = sprintf('%s_intFrequ_%s_modScore_%s',   img_outfile, num2str(round(freqPeak,2)) , num2str(round(mod_score,2)));
     print(fig1, '-dpdf', strcat(img_outfile,'.pdf'),'-painters');
     close
     end
     
     
end

function [cor, lag, smooth] = SmoothCor(cor, lag, t_bin)
    
    std_smooth_kernel = .005;
    
    kernel = pdf('normal', -std_smooth_kernel*10/t_bin:std_smooth_kernel*10/t_bin, 0, std_smooth_kernel/t_bin); % convolve with gaussian

    if isempty(cor), smooth = []; return; % if no signals, then send it back
    
    else
        
        smooth = zeros(size(cor));
        
        for i = 1:size(cor,2)  

            smooth(:,i) = ndnanfilter(cor(:,i), kernel(:), []);
        
        end
        
        lag = lag(round(end/2):end,:);
        
        cor = cor(round(end/2):end,:);
        
        smooth = smooth(round(end/2):end, :);

    end 
    
end

function [S, f] = getChronuxSpectrum(cor, t_bin)

    params.tapers = [1.25, 1];
    params.fpass = [];
    params.Fs = 1/t_bin;
    params.pad = 6;  

    [S,f]=mtspectrumc(cor,params); % compute spectrum of cor

end