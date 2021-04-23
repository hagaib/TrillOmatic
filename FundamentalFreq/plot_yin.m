function h = plot_yin(f0 , dips , time , axes , colors)
% Plots yin f0 results according to probabilities
% best : below 0.2% periodicity
% good : 0.2% < good < 0.4%
% bad > 0.4%

% f0 , dips , time : output of yin pitch estimation algorithm.
% axes: where to draw
% colors: optional. in the form [c1 , c2 , c3]. c1= best color , c2= good
% color , c3 = bad color


if nargin < 5
    colors = ['g' , 'c' , 'm'];
end
hold (axes , 'on');
best = dips < 0.2;
time_plot = time;
time_plot(~best) = nan;
f0_plot = f0;
f0_plot(~best) = nan;
p1 = plot(time_plot , f0_plot , sprintf([colors(1) , '.-']) , 'LineWidth' , 3);

good = (dips < 0.4) & (dips >= 0.2);
time_plot = time;
time_plot(~good) = nan;
f0_plot = f0;
f0_plot(~good) = nan;
p2 = plot(time_plot , f0_plot , sprintf([colors(2) , '.-']), 'LineWidth' , 3);


time_plot = time;
% time_plot(good | best) = nan;
f0_plot = f0;
f0_plot(good | best) = nan;
p3 = plot(time_plot , f0_plot , sprintf([colors(3) , '.-']), 'LineWidth' , 3);

h = [p1 ;p2 ;p3];
hl = legend(h ,'best','good' , 'else');
h = [h ; hl];

legend('show')

hold (axes , 'off');