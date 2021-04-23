function log = segments2logical(segments  , time , fs)
%returns logical array of length "length(time)" corresponding to time
%information in segments matrix;
% segments is a matrix with dim(segments,1)=2

log = zeros(size(time));

if(size(segments,1)~=2)
    segments = segments';
    if(size(segments,1)~=2)
        error('dim(segments,1)~=2 for function to work correctly');
    end
end

segments = floor(1+(segments-time(1))*fs);
for i=1:size(segments,2)
    log(segments(1,i):segments(2,i)) = 1;
end

log = logical(log);
