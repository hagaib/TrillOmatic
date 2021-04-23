function [meanprob , badprecent] = mean_yinprobs_trill(yin , min_segment_time , periodic_bw)
% approximates mean yin dip and precentage of high dips (over 0.4).
% precentage is not calculated out of entire signal, but only a portion of
% the signal which is believed to be highly periodic

if(nargin<2)
    min_segment_time = 0.018;
    periodic_bw = [2000 , 3500];
end
[periodictime , nonperiodictime , persegs , nonpersegs] = pertime(yin.dips , yin.time , min_segment_time ,yin.f0, periodic_bw);
itime = longtrill_bulk_of_mass(persegs , yin , 1.5 , 0);

probs = yin.dips(itime(1):itime(2));
probs(isnan(probs)) = 1;
meanprob = mean(probs);
badprecent = sum(probs(probs>0.4))/length(probs);