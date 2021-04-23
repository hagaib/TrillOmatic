function [lbout , hbout] = clip_to_min_bw( lowb , highb , minBW )
%CLIP_TO_MIN_BW Summary of this function goes here
% higherb = higher bound of BW
% lowerb = lower bound of BW
% minBW = minimum bandwidth allowed

if(highb-lowb) < minBW
    centerf = 0.5*(highb+lowb);
    hbout = centerf + minBW/2;
    lbout = centerf - minBW/2;
else
    lbout = lowb;
    hbout = highb;
end
