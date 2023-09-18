function [t, x, y, angle,n_samples] = ProcessVideoData(file,diode)
%diode: 0 - luminance, 1 - doublediode (processed by a separate function),
%2 - red, 3 - green

handles = readVideoData(file); %displays number of records
[dTargets,trackingColour] = decodeTargets(handles.targets);

if ~exist('diode','var')||isempty(diode)
    diode = 0;
end
switch diode
    case 1 % double diode, process with a separate function
        [t, x, y, angle, n_samples] = ProcessVideoData_double(file);
        return
    case 2 % Use red signal
        [x,y,targets,exFlag] = extractPosition(dTargets,[0 1000],[0 1000],2);
        x = x'; y = y'; t = handles.post;color = ('red');
    case 3 % Use Green Signal
        [x,y,targets,exFlag] = extractPosition(dTargets,[0 1000],[0 1000],3);
        x = x'; y = y'; t = handles.post;color = ('green');
    otherwise % Use luminance
        [x,y,targets,exFlag] = extractPosition(dTargets,[0 1000],[0 1000],1);
        x = x'; y = y'; t = handles.post;color = ('lum');
end

%% If videoFix screwed up timestamps, (makes 10^10 instead 10^4)
if t(1) < 0
errordlg(sprintf('Negative timestamps in: %s', file));
end
if t(1) > 10^14
   errordlg('Fixing Video timestamps')
   t = t * 10^-6;
end   

n_samples = length(x);
%% Suppress postions at [0,0] because this is caused by bad tracking
%[x,y,t] = suppressZeros(x,y,t);
    
ind = x==0;
%t(ind) = NaN; %[];
x(ind) = NaN; %[];
y(ind) = NaN; %[];

%if sum(~(isnan(x)))<100; errordlg(fprintf('Insufficient videodata, please check')); end

%% Speedfilter, too fast is trackingerror:  %introduced @MH 2/2/2018
filter_rate = 100; %this is roughly in cm / sec based on current camera setting of 2pix roughly 1cm
[x,y] = speedfilter_MH(x,y,t,filter_rate);
%plot(rates,'r');hold on

%% Remove reflections  %introduced @MH 02/03/2018
% if successfull no speed threshholding needed anymore
% find_outside_pixels(x, file)  %check if there are still values outside
% find_outside_pixels(y, file)  %check if there are still values outside
 


sampRate = 30;
timeThreshold = 0.5;
 [x,y] = interpolatePosition_sl(x,y,timeThreshold,sampRate)
 
%% Fix NaNs
%x = fix_nans(x);
%y = fix_nans(y); 
%{
figure, plot(x,y)
title('Before reflections removal')

x(x<frame_begin) = NaN;
x(x>frame_end) = NaN;
[frame_begin, frame_end] = find_outside_pixels(y);
y(y<frame_begin) = NaN;
y(y>frame_end) = NaN;
figure, plot(x,y)
title('Removed reflections')
%}


[x,y] = meanpath(x,y);
figure
plot(x,y,'g');
title('After speedfiltering jumps and fixing NaNs');
%findpeaks(rates,30,'Annotate')
%test_x = diff(x);
%test_y = diff(y);
%test_xy = test_x.*test_y;
%plot(test_x,'b');hold on
%plot(test_y,'g');hold on

%% Smoothening of positioning?
%[x,y] = interpolatePosition_sl(x,y,timeThreshold,sampRate);
          
%% Add angle
angle(1:15) = NaN;
for i = 16:(length(x)-15)
    %disp(i);
    if (x(i+15) == x(i-15)) && (y(i+15) == y(i-15))
        angle = NaN;
    else
        [a, d] = cart2pol(x(i+15)-x(i-15),y(i+15)-y(i-15));
        angle(i) = a;
    end
end 
angle((length(x)-14):length(x)) = NaN;

%{
%% Make figure
fh = figure();
switch diode
    case 1
    plot(x,y,'b');     
    case 2    
    plot(x,y,'r'); 
    case 3
    plot(x,y,'g');
    otherwise
    plot(x,y,'k');
end  

sequ = 1; fname = fullfile(file(1:end-7), sprintf('runVelocity_Tracking_%d.bmp', diode));

%%Use this in case you want several sequential files in a folder
%{
while exist(fname, 'file')
sequ = sequ+1;fname = fullfile(file(1:end-7), sprintf('runVelocity_Tracking_%d.bmp', sequ));
end
%}

saveas(fh,fname);
saveas(fh,fullfile(file(1:end-7), sprintf('runVelocity_Tracking_%d.ai', diode)),'epsc');
%hold on
 close(fh);
 %}

% end
%figure(2); for i=1:100:length(gx)-100 hold off; plot(gx(i:i+100),gy(i:i+100),'g.'); hold on; plot(rx(i:i+100),ry(i:i+100),'r.'); plot(posx(i:i+100),posy(i:i+100),'b.'); axis([300 700 0 400]); pause; end

function [x,y] = meanpath(x,y)

temp_x = x;
temp_y = y;
for cc=2:length(x)-2
  x_window = x(cc-1:cc+1); y_window = y(cc-1:cc+1);
  temp_x(cc) = nanmean(x_window); temp_y(cc) = nanmean(y_window);
end

x = temp_x;
y = temp_y;

function [x,y] = interpolatePosition_sl(x,y,timeThreshold,sampRate)

sampThreshold = floor(timeThreshold * sampRate);

index_isvalid = find(~isnan(x) & ~isnan(y));
count_NaNs = 0;
if numel(index_isvalid)<2
    errordlg('Position data empty! Please check if splitting data was successfull');
end
for i = index_isvalid(1):index_isvalid(end)
    if isnan(x(i)) || isnan(y(i)) 
        count_NaNs = count_NaNs+1;
    else
        if count_NaNs > 0 && count_NaNs < sampThreshold
            for j = last_valid+1:i-1
                x(j) = (x(i)*(j-last_valid)+x(last_valid)*(i-j))/(i-last_valid);
                y(j) = (y(i)*(j-last_valid)+y(last_valid)*(i-j))/(i-last_valid);
            end
        end
    last_valid = i;
    count_NaNs = 0;
    end
end

function [ ] = find_outside_pixels(x,file)
all_values = unique(x(:));
stepsize = diff(all_values);
if max(stepsize)>2
    %errordlg(sprintf('Error processing\n %s\n Please videofix',file));
fileID = fopen('C:\Users\Silvia\Dropbox\SilviaProjectCode\ToDo_Videofix.txt','a');
fprintf(fileID,'\n%s\n', file);
fclose(fileID);
end
%{
continuous_values = stepsize<2;
discontinuation = find(~(continuous_values));
discontinuation = [1; discontinuation; numel(all_values)];  %add beginning and end
[biggest, biggest_init_idx] = max(diff(discontinuation));
frame_begin = all_values(discontinuation(biggest_init_idx)+1);
frame_end = all_values(discontinuation(biggest_init_idx+1));
%}


function [x,y] = speedfilter_MH(x,y,t,filter_rate)
a = hypot(0.5*diff(x),0.5*diff(y));
pertime = diff(t)*10^-6;  %conversion to seconds 
rates = a./pertime;
filter_values = find(rates>filter_rate);
%% Use this to plot
figure
scatter(x(filter_values),y(filter_values),'r*');hold on
plot(x,y,'b'); hold off
title('Before speedfiltering jumps');
%% Filter now points marked in red
x(filter_values) = NaN;
y(filter_values) = NaN;

function [x] = fix_nans(x)
fprintf('Fixing missing values in x and y\n')
% Fix these:
to_fix = find(isnan(x));
while ~isempty(to_fix)  % check if there are still NaNs
do_now = to_fix(diff(to_fix)>1); %fix the easy ones first, the others have to work in sort of iterations
        if ((numel(do_now))<2)
            do_now = [do_now,to_fix(1), to_fix(end)];
        end    
for fff = do_now   %% these are the values to fix
   fprintf('Fixing: %s\n', num2str(fff));
  % do parametersweep till interpolation can be done
  searchradius = 0;
  while isnan(x(fff))
  searchradius = searchradius + 1;    
  bound1 = max([(fff-searchradius), 1]); % make sure not to go below 0
  bound2 = min([(fff+searchradius), numel(x)]); % make sure not to go beyond x(end)
  fprintf('Interpolating from %s to %s\n', num2str(bound1),num2str(bound2));
  if (~isnan(bound1) && ~isnan(bound2))  %check that none of the two is a NaN value
  x(fff) = nanmean(x(bound1:bound2)); % interpolate the value
  elseif bound1==1 && ~isnan(bound2)
  x(fff) = nanmean(x(bound1:bound2)); % interpolate the value
  elseif ~isnan(bound) && bound2==numel(x)
  x(fff) = nanmean(x(bound1:bound2)); % interpolate the value    
  end
  if searchradius > 100 %if it becomes ridiculous stop it
    break
  end
  
  end
end

to_fix = find(isnan(x)); % check if there are still NaNs
end