function yy = prcfilt(xx , winsize , prc , nanflag)
%Precentile Filter - same as median filter but result is the "prc"
%precentile in each window. (if prc==50 , result is same as median filter)
%winsize must be an integer (number of samples)
% nanflag: 'omitnan' || 'includenan'
% data is zero padded
tic
if(nargin < 4)
    nanflag = 'includenan';
end

xx = xx(:);
L = length(xx);
transient = floor(winsize/2);

yy = nan(size(xx));

xx = [zeros(transient ,1) ; xx ; zeros(transient+1 ,1)];

nanarray = isnan(xx);
if(isequal(nanflag ,'includenan'))
    nanarray = zeros(size(xx));
end

countnan = sum(nanarray(1:winsize-1));

for i=1:L
    winvals = xx(i:i+winsize-1);
    yy(i) = percentile(winvals , prc , countnan);
    countnan = countnan - nanarray(i) + nanarray(i+winsize);
end
toc

% Percentile calculated the k'th percentile of x. This function is similar 
% to, but generally much faster than MATLAB's prctile function.
function y = percentile(x, k , countnan)
x = sort(x);
n = size(x,1) - countnan;

p = 1 + (n-1) * k / 100;

if p == fix(p)
    y = x(p);
else
    r1 = floor(p); r2 = r1+1;
    y = x(r1) + (x(r2)-x(r1)) * k / 100;
end