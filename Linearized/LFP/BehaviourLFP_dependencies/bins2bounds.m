function [bounds] = bins2bounds(bins)
%Defines boundaries of areas / binnning
%   Detailed explanation goes here

bounds = [0, bins(1), sum(bins(1:2)), sum(bins(1:3)), sum(bins(1:4)), sum(bins(1:5))];  

end

