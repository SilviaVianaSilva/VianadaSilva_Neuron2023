function [] = plot_Behav_quietSl_MovSl(filename,f,behav,sleepscored,quiet,active);        
%To check if z-scoring to quiet periods in S1 in Box works well for every
%animal
% Plot the average power during behaviour on the box, during the quiet period and
% active period in the box; sleep z-scored will be on different scale, so
% excluded for now.
% Called by ConvertRegionalEEGs
%
%   Matthias Haberl, Nov 3rd, 2017

figure()
b = plot(f, nanmean(behav,2),'.-k'); hold on
q = plot(f, quiet,'.-g'); hold on
a = plot(f, active,'.-r'); hold off

legend([b,q,a],{'behaviour','quiet box','active box'})
set(gca, 'xscale', 'log', 'XTick', [2, 4, 8, 30, 100,300], 'Box', 'off', 'FontName', 'Arial', 'FontSize', 12, 'TickDir', 'out');
ylabel('Mean Power');xlabel('Frequency');
print('-dpdf', filename,'-opengl') %-painters :painters is vectorgraphics but slow and too big
saveas(gcf,filename,'epsc')

figure()
s = plot(f, mean(sleepscored,2),'.-g');
legend([s],{'sleepscored'})
set(gca, 'xscale', 'log', 'XTick', [2, 4, 8, 30, 100,300], 'Box', 'off', 'FontName', 'Arial', 'FontSize', 12, 'TickDir', 'out');
ylabel('Power z-scored to Rest 1');xlabel('Frequency');
title('Power z-scored to resting phase 1');
filename2 = [filename,'sleepscored'];
print('-dpdf', filename2,'-opengl') %-painters :painters is vectorgraphics but slow and too big
saveas(gcf,filename2,'epsc')
close all

end

