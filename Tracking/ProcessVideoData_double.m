function [t, x, y, angle, n_samples] = ProcessVideoData_double(file,outdir);

%% Load data within Range
borderX = [0 1000];
borderY = [0 1000];
sampRate = 30;
timeThreshold = .5;
%speedlimit = 100; %in pixels/sec 

handles = readVideoData(file) %displays number of records
[dTargets,trackingColour] = writeVelocity.decodeTargets(handles.targets);

[rx,ry,rtargets,exFlag] = extractPosition(dTargets,[0 1000],[0 1000],2);
[gx,gy,gtargets,exFlag] = extractPosition(dTargets,[0 1000],[0 1000],3);
t = handles.post; %t = t';

%% If videoFix screwed up timestamps, (makes 10^10 instead 10^4)
if t(1) < 0
errordlg(sprintf('Negative timestamps in: %s', file));
end
if t(1) > 10^14
   errordlg('Fixing Video timestamps')
   t = t * 10^-6;
end

n_samples = length(rx);


%% Remove zeros
ind = find(rx == 0);
rx(ind) = NaN; %[];
ry(ind) = NaN; %[];
ind = find(gx == 0);
gx(ind) = NaN; %[];
gy(ind) = NaN; %[];


%% Speedfilter, too fast is trackingerror:  %introduced @MH 2/2/2018
filter_rate = 100; %this is roughly in cm / sec based on current camera setting of 2pix roughly 1cm
[gx,gy] = speedfilter_MH(gx,gy,t',filter_rate);
[rx,ry] = speedfilter_MH(rx,ry,t',filter_rate);
%plot(rates,'r');hold on

[rx,ry] = interpolatePosition_sl(rx,ry,timeThreshold,sampRate);
[gx,gy] = interpolatePosition_sl(gx,gy,timeThreshold,sampRate);

%if red and greed or more than 2 std apart than the mean, one of the
%tracking points is wrong; exclude from calculating xy coordinate
grdist = sqrt((rx-gx).*(rx-gx)+(ry-gy).*(ry-gy));
gr_factor = 2; %number of std
ind = find(grdist > nanmean(grdist)+nanstd(grdist)*gr_factor);

rx(ind) = NaN; %[];
ry(ind) = NaN; %[];
gx(ind) = NaN; %[];
gy(ind) = NaN; %[];

[rx,ry] = interpolatePosition_sl(rx,ry,timeThreshold,sampRate);
[gx,gy] = interpolatePosition_sl(gx,gy,timeThreshold,sampRate);

%% Speedlimit
%[rx,ry,rtargets,exFlag] = writeVelocity.removeJumps_sl(rtargets,[0 1000],[0 1000],2,speedlimit);
%[gx,gy,gtargets,exFlag] = writeVelocity.removeJumps_sl(gtargets,[0 1000],[0 1000],3,speedlimit);


%% Remove reflections  %introduced @MH 02/03/2018
% if successfull no speed threshholding needed anymore %for now used to
% indicate if there is additional fix to do
% find_outside_pixels(gx,gy,rx,ry, file);  %check if there are still values outside
 
%% Fix NaNs  %use this or interpolation a la sl
%rx = fix_nans(rx);
%ry = fix_nans(ry);
%gx = fix_nans(gx);
%gy = fix_nans(gy);

%% Smoothening
[gx,gy] = meanpath(gx,gy);
[rx,ry] = meanpath(rx,ry);

%
rx = double(rx); ry = double(ry);
gx = double(gx); gy = double(gy);

%% Merge colors to single track
x = (rx+gx)/2; x = x';
y = (ry+gy)/2; y = y';


fh = figure(); plot(gx,gy,'g'); hold on; plot(rx,ry,'r'); plot(x,y,'b'); 
save_dir = fullfile(outdir,'TrackingImages');
mkdir(save_dir);
saveas(fh,strcat(save_dir,'Interpolated_Tracking_before_rescaling.png')); 
close(fh);
%figure(2); for i=1:100:length(gx)-100 hold off; plot(gx(i:i+100),gy(i:i+100),'g.'); hold on; plot(rx(i:i+100),ry(i:i+100),'r.'); plot(posx(i:i+100),posy(i:i+100),'b.'); axis([300 700 0 400]); pause; end

%if red and greed are more or less than 1.5 std from the mean, angle is likely 
%inaccurate; exclude from calculating angle
ind = find(grdist > nanmean(grdist)+nanstd(grdist)*gr_factor | grdist < nanmean(grdist)-nanstd(grdist)*gr_factor);
rx(ind) = NaN; %[];
ry(ind) = NaN; %[];
gx(ind) = NaN; %[];
gy(ind) = NaN; %[];

angle = cart2pol(rx-gx,ry-gy); angle = angle'; 


function [x,y] = interpolatePosition_sl(x,y,timeThreshold,sampRate)

sampThreshold = floor(timeThreshold * sampRate); % do not fill NaNs if more then this (30*0.5  =15)

index_isvalid = find(~isnan(x) & ~isnan(y));
count_NaNs = 0;
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

            
    


function [x,y] = meanpath(x,y)

temp_x = x;
temp_y = y;
for cc=2:length(x)-2
  x_window = x(cc-1:cc+1); y_window = y(cc-1:cc+1);
  temp_x(cc) = nanmean(x_window); temp_y(cc) = nanmean(y_window);
end

x = temp_x;
y = temp_y;

function [ ] =  find_outside_pixels(gx,gy,rx,ry, file)
all_values_x = unique([gx(:), rx(:)]);
stepsize_x = diff(all_values_x);
all_values_y = unique([gy(:), ry(:)]);
stepsize_y = diff(all_values_y);

if (max(stepsize_x)>2) || (max(stepsize_y)>2)
%    errordlg(sprintf('Error processing\n %s\n Please videofix',file));
%% to print to a text file if jumps above 2
fileID = fopen('C:\Users\Silvia\Dropbox\SilviaProjectCode\ToDo_Videofix_doublediode.txt','a');
fprintf(fileID,'%s\r\n', file);
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
close
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
  end
  if searchradius > 100 %if it becomes ridiculous stop it
    break
  end
  
  end
end

to_fix = find(isnan(x)); % check if there are still NaNs
end