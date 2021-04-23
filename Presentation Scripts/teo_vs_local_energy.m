filenames =  get_list_csv('filenames.csv');

filename = filenames{9};


[x , fs] = audioread(filename);
x = x(:,1);
N= length(x);
t = 1:N;
t = (t-1)/fs;



[b,a] = butter(7, [1800 , 3500]/fs*2);
xfilt = filtfilt(b,a,x);

energy = hamming_energy(xfilt , fs , 0.002 , 1);

teo = teager(xfilt);
teo = [teo(1) ; teo ; teo(end)];

L = floor(0.01*fs); 
if(~mod(L,2)) , L=L+1; end
transient = (L-1)/2;
g = gausswin(L);
g = g/sum(g);

teofilt = filter(g , 1 , teo);
teofilt = [teofilt(transient+1:end) ; zeros(transient , 1)];

E_g = filter(g , 1 , energy);
E_g = [E_g(transient+1:end) ; zeros(transient , 1)];


%% PRI Estimation

fourierteo = fft(teofilt);
fourierenergy = abs(fft(energy));

fourierabs = abs(fourierteo);
[m1 , f1]  = max(abs(fourierabs));
[m2, f2] = max(abs(fourierenergy));

[p , tpeaks] = findpeaks(fourierabs(1:floor(N/2)));
[~ , maxt] = max(p);
maxt = tpeaks(maxt);
PRI = 1/((maxt-1)/N*fs)



theta = angle(fourierteo(maxt));
amp = abs(fourierteo(maxt))/N*6;


%% Drawing Everything
figure(3)
p1 = subplot(3,1,1);
plot(t,x);
title(filename)
legend('Original Signal')

p2 = subplot(3,1,2);
plot( t , teofilt, 'g')
% plot(t, energy , 'b' , t,teo ,'r'  , t , teofilt, 'g');
hold on
plot(t , mean(teofilt) + amp*cos(t*2*pi*fs*(maxt-1)/N+ theta) , 'k')
hold off
legend('Smooth Energy' , 'Sinusoidal Approximation');
xlabel('Time')


p3 = subplot(3,1,3);
plot((0:(N-1))/N*fs , abs(fourierteo) )%, (0:(N-1))/N*fs,abs(fourierenergy));
legend('Fourier TEO' , 'Fourier Energy')
xlabel('Hertz')

linkaxes([p1,p2],'x')

figure(4)
plot(t,xfilt,'b')
hold on
plot(t, energy , 'b' , t,teo ,'r'  , t , teofilt, 'g');
hold off

figure(5)
plot(t,xfilt,'b')
hold on
plot(t, 20*energy , 'k' ,'Linewidth' , 2)
plot(t , 20*E_g , 'g' , 'Linewidth' , 2)
hold off
legend('Bandpass Filtered Signal' ,'Energy' ,  'Gaussian Smoothed Energy')
title(filename)