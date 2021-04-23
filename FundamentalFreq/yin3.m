function [ f0 , dips , time ] = yin3(signal , fs , windur , threshold , minfreq, stepsize)
%Approximates fundamental frequency of signal according to the celebreated
%YIN algorithm
Tmax = floor(1/minfreq * fs) + 2;
winsize = floor(windur * fs);
hops = winsize:stepsize:length(signal)-winsize-Tmax+1;

low_pass = fir1(48 , 0.5);
signal = filtfilt(low_pass , 1 , signal);
dd = diff_function(signal , Tmax , winsize , hops);
[dd , ddx] = parabolic_interpolation(dd);
dd = cummnorm(dd);


f0 = nan(size(hops));
dips = nan(size(hops));
periods = nan(size(hops));

for i = 1:length(hops)
    hop = hops(i);
    d = dd(hop,:);
    per = period_from_diff(d ,[1 Tmax] , threshold , 2);
%     if(hop/fs>5.58)
%         ttt = hop/fs
%         figure(2); plot(d);
%         fff = fs/per
%     end
    if ~isnan(per)
        dips(i) = d(per);
        periods(i) = per;
        f0(i) = fs/(per+ddx(hop, per));
    end
    
end


range = -floor(0.5*Tmax):ceil(0.5*Tmax);
temp_periods = nan(size(periods));
for i=1:length(hops)
    interval = i + range;
    interval = interval(interval>0);
    interval = interval(interval <= length(hops));
    [~ , j] = min(dips(interval));
    j = interval(1) +j-1; 
    temp_periods(i) = periods(j);
end
periods = temp_periods;



low = 0.6;
high= 1.5;
for i = 1:length(hops)
    if(i==857)
    end
    hop = hops(i);
    per = periods(i);
    dip = dips(i);
%     if i >= 1100
%         per
%         i
%     end
    if (~isnan(per))
%         for j=hop+range%shop around the interval
%             if j == hop , continue; end
        d = dd(hop,:);
        freq_range = [floor(low*per) , ceil(high*per)];
%         if(min(d)<1)
%             figure(1)
%             plot(d)
%             hop
%         end
        refined_per = period_from_diff(d , freq_range , 0 , 1);
        if ~isnan(refined_per)
            dip = d(refined_per);
            per = refined_per + ddx(hop, refined_per);
            periods(i) = per;
%             if refined_dip < dip
%                 per = refined_per;
%                 dip = refined_dip;
%                 best_time = j;
%             end
%         end
        else
            periods(i) = nan;%per + ddx(hop , per);
            dip = nan;
        end
        dips(i) = dip;
    end    
end

% f0_refined = fs./periods;
f0 = fs./periods;
time = hops / fs;

function per = period_from_diff(d , range , threshold , ver)
dip = inf;
per = nan;
if ver == 2 , threshold = threshold+min(d); end
range = range(1):range(2);
if range(1) == 1 , range = range(2:end); end
if (range(end) >= length(d)) , range = range(1):length(d)-1; end
for i=range
    if(d(i) < d(i-1) && d(i) < d(i+1))
        if d(i) < threshold
            per = i;
            return;
        elseif d(i) < dip
            per = i;
            dip = d(i);
        end
    end
end


function cumdd = cummnorm(dd)
cumdd = ones(size(dd));
cumsum = dd(:,1);

for tau=2:size(dd,2)
    cumsum = cumsum + dd(:,tau);
    cumdd(:,tau) = dd(:,tau) ./ cumsum;
    cumdd(:,tau) = cumdd(:,tau) * tau;
end

%cumulative mean-normalize 
function y=cumnorm_de_chev(x)
x=x;
[m,n]=size(x);
y = cumsum(x);
y = (y)./ (eps+repmat((1:m)',1,n)); % cumulative mean
y = (eps+x) ./  (eps+y);


function dd = diff_function(signal , Tmax , winsize , hops)
if nargin < 3 , hops = nan; end
dd = zeros(length(signal) , Tmax);
for tau=1:Tmax
    sig_diff_sq = signal(1:end-tau) - signal(1+tau:end);
    sig_diff_sq = sig_diff_sq .* sig_diff_sq;
    dd(1, tau) = sum(sig_diff_sq(1:winsize));
    if isnan(hops)
        for t=2:length(sig_diff_sq)-winsize+1
            dd(t, tau) = sum(sig_diff_sq(t:t+winsize-1));
        end
    else
        for t=hops
            dd(t, tau) = sum(sig_diff_sq(t:t+winsize-1));
        end
    end
end

function [dd , ddx] = parabolic_interpolation(dd)
ddx = zeros(size(dd));
for j=1:size(dd,1)
    d = dd(j,:);
    dx = zeros(length(d),1);
    to_interpolate = nan(size(dx));
    for i=2:length(d)-1
        if(d(i)<d(i-1) && d(i)<d(i+1))
            to_interpolate(i) = i;
        end
    end
    to_interpolate = to_interpolate(~isnan(to_interpolate));
    for i=1:length(to_interpolate)
        index = to_interpolate(i);
        [dx(index),d(index)] = parabolic_min(d(index-1:index+1));
    end
    dd(j,:) = d;
    ddx(j,:) = dx;
end


function [min_x , min_y] = parabolic_min(y)
    if (length(y) ~= 3)
        error('Must provide 3 points for parabolic interpolation');
    end
    a = 0.5*(y(1) + y(3)) - y(2);
    b = 0.5*(y(3) - y(1));
    c = y(2);
    
    min_x = - b/a/2; %x_min = -b/2a
    min_y = c - b*b/a/4; %y_min = c-b^2/4a

