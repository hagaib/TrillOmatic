function y = prctile(v, p)
% percentile p of data points v
pctl = @(v,p) interp1(linspace(0.5/length(v), 1-0.5/length(v), length(v))', sort(v), p*0.01, 'spline');
y = pctl(v, p);