function y = moving_average(x,n , start , stop)
%returns non periodic moving average of signal 
%averages around sample x(i) in the domain [i-n,i+n]
%   y = averaged signal
%   x = original signal
%   n = half the "width" of the average
%   start, stop : optional parameters. in case only a portion of original
%   signal needs to be averaged
%   sample = optional parameter. average only around specific sample
    
    if (nargin < 3)
        start = 1;
        stop = length(x);
    end 
    
    trans = 0;
    if (2*n+1 > length(x)) , error('moving_average:not_enough_samples' , 'signal is too short for desired average'), end
    if (size(x,2) == 1) , x = transpose(x); trans=1; end
    if (size(x,1) ~= 1) , error('input must be 1 dimensional signal'), end
    if (stop < 1 || start < 1 || stop > length(x) || start > length(x) || stop < start) 
        error('illegal start or stop values'), end
    
    y = zeros(size(x));
    
    i=1;
    if i < start
        y(i:start-1) = x(i:start-1);
    end
    
    y(start) = sum(x(max(1,start-n):start+n));
    i = start+1;
    while i<=n+1
        y(i) = y(i-1) + x(i+n);
        i=i+1;
    end
    
    while i <= min(length(x) - n , stop)
        y(i) = y(i-1) + x(i+n) - x(i-n-1);
        i=i+1;
    end
    
    while i <= stop
        y(i) = y(i-1) -x(i-n-1);
        i=i+1;
    end
    
    if i <= length(x)
        y(i:end) = x(i:end);
    end
    
    indx = start:n;
    if (~isempty(indx)) , y(indx) = y(indx) ./ (indx+n); end
    indx = length(y)-n+1:stop;
    if (~isempty(indx)) ,  y(indx) = y(indx) ./ flip(2*n-length(indx)+1:2*n); end
    indx = max(start,n+1):min(stop,length(y)-n);
    if (~isempty(indx)) , y(indx) = y(indx) / (2*n+1); end
    
    if(trans) , y= transpose(y); end
end