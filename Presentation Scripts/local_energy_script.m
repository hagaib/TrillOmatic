M = 10001;
t = linspace(-20*pi , 20*pi , M);
y = (cos(t/6)) .* cos(t);
L = 1501;
w = hanning(L);

support = (M-1)/2+(-(L-1)/2:(L-1)/2);

ywin = zeros(size(y));
ywin(support) = y(support) .* w';

enwin = conv(y.^2 , w.^2);

figure(2)
p1 = subplot(3,1,1);
plot(t,y);
hold on
plot(t(support) , w);
hold off
legend('cos(t/6)\cdotcos(t)' , 'Hanning Window')

p2 = subplot(3,1,2);
plot(t , ywin)
legend('Windowed Function for t=0')

p3 = subplot(3,1,3);
plot(t(1+(L+1)/2:end-(L+1)/2) , enwin(L+1:end-L));
legend('Local Energy for all t in range')

linkaxes([p1,p2,p3],'x')