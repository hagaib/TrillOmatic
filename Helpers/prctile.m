function y = prctile(x, k)
% percentile p of data points v
x = sort(x);
n = size(x,1);

p = 1 + (n-1) * k / 100;

if p == fix(p)
    y = x(p);
else
    r1 = floor(p); r2 = r1+1;
    y = x(r1) + (x(r2)-x(r1)) * k / 100;
end