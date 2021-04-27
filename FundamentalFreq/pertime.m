function [periodictime , nonpertime , persegs , nonpersegs] = pertime(dips , time , mindur ,  f0 , periodic_bw, perprob , energyvad)
% calculates segments of signal in which it is highly periodic , and times
% in which it is highly non-periodic. 
% calculation is based on yin output.
% energyvad = vad energy binary detection vector from prior vad (optional)
%f0 ,periodic_bw : optional parameters . In case only frequencies of
%specific bandwidh are relevant

% output:
% periodictime: vector of indicator function for periodic portions of
% signal
% nonpertime: vector of indicator function for non periodic portions of
% signal
% persegs: table of starting and ending times for periodic portions of
% signal
% nonpersegs: table of starting and ending times for non periodic portions
% of signal.


if(nargin < 6) 
    perprob = 0.2;
end

periodictime = dips < perprob;

if(nargin > 3)
    periodictime = periodictime & f0>periodic_bw(1) & f0<periodic_bw(2); 
end

if(nargin > 6) 
    periodictime = periodictime & energyvad;
end

periodictime = logical(medfilt1(double(periodictime),3));
periodictime(1) = 0;
periodictime(end) = 0;
dpertime = periodictime(2:end) - periodictime(1:end-1);
startsegs = find(dpertime>0);
endsegs = find(dpertime<0);

segdurs = time(endsegs) - time(startsegs);

mindur = max(mindur , prctile(segdurs , 50))
% mindur = prctile(segdurs , 70)
segdurs = segdurs < mindur; % ignore segments shorter than mindur

persegs = [time(startsegs(~segdurs)) ; time(endsegs(~segdurs))];
startsegs = startsegs(segdurs);
endsegs = endsegs(segdurs);
for i=1:length(startsegs)
    periodictime(startsegs(i):endsegs(i)) = 0;
end

nonpertime = dips > 0.4;
nonpertime = medfilt1(double(nonpertime) , 3);
nonpertime(1) = 0;
nonpertime(end) = 0;
dnonpertime = nonpertime(2:end) - nonpertime(1:end-1);
startsegs = find(dnonpertime>0);
endsegs = find(dnonpertime<0);
segdurs = time(endsegs) - time(startsegs);
segdurs = segdurs < mindur/2;

nonpersegs = [time(startsegs(~segdurs)) ; time(endsegs(~segdurs))];
startsegs = startsegs(segdurs);
endsegs = endsegs(segdurs);
for i=1:length(startsegs)
    nonpertime(startsegs(i):endsegs(i)) = 0;
end
end