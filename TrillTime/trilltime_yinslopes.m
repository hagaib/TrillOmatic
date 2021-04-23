function [time , yinslopvar] = trilltime_yinslopes(yin , windur , threshold , restricttime)
% find global trill time according to yinslopevar filter
%   time = trilltime_yinslopes(yin , wint , 100 , [0.56 , 1.74]);


yinslopvar = yinslopevarfilt(yin , windur);
tt = yin.time(sqrt(yinslopvar) < threshold);
tt = tt(tt>restricttime(1) & tt<restricttime(2));
time = [tt(1) , tt(end)] + 0.5*windur*[-1,1];

