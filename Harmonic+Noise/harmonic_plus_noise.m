function yy = harmonic_plus_noise(xx , fs , freq ,windur , stepdur)

yy = zeros(size(xx));

winsize = floor(windur*fs);
if mod(winsize , 2) , winsize = winsize+1; end 
hann = hanning(winsize);
stepsize = floor(stepdur*fs);
L = 4;


for i=winsize/2:stepsize:(length(xx)-winsize/2)
    samples = i-winsize/2+1:i+winsize/2;
    win = xx(samples);
    
%     f = abs(fft(hann.*win)/winsize);
%     P = f(1:winsize/2+1);
%     P(2:end-1) = P(2:end-1) + wrev(f(winsize/2+2:end));
%     plot((0:winsize/2)*fs/winsize,P);
    
    B = repmat(samples'/fs,1, 2*L+1);
    B = B * freq;
    B = B * diag(-L:L);
    B = exp(2*pi*1i*B);
    
    R = (B)'*diag(hann.^2)*B;
    b = (B)'* diag(hann.^2)*win;
    x = R\b;
    yy(samples) = hann.*real(B*x)+yy(samples);
end