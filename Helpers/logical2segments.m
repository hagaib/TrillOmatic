function segments = logical2segments( logical , fs )
%LOGICAL2SEGMENTS returns a 2xsegcount matrix of segments (beginning and
%end) in seconds or samples
%logical: is a vector of zeros and ones correspoding to desired segments 
%fs (optional): is sample rate. if it's not specificed , output will be in
%samples. otherwise, it is in seconds

logical = logical(:)';

dl = logical(2:end) - logical(1:end-1);
start = 1 + find(dl == 1);
stop = find(dl == -1);

if(start(1) > stop(1)) 
    start = [1 start];
end
if(stop(end) < start(end))
    stop = [stop length(logical)];
end

segments = [start ; stop];
if(nargin==2)
    segments = (segments-1)/fs;
end