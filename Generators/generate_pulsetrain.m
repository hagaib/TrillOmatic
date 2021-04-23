function yy  = generate_pulsetrain(pulses , fs , PRI , SNR)
%Generates pulse train.
% pulses = cell array
% PRI = every pulse is generated at PRI interaval
%  SNR is defined as 20*log10(Esignal/Enoise)
if (nargin < 4) , SNR = inf;  end 


yy = zeros(1 , floor(fs*PRI)*length(pulses));

total_length = 0;
for i=1:length(pulses)
    pulse = pulses{i};
    
    start = floor(fs*PRI * i);
    
    window = hanning(length(pulse))';
    
    total_length = length(pulse) + total_length;
    
    
    yy(start:(start+length(pulse)-1)) = pulse .* window;

end

if(SNR < inf)
    noise = randn(total_length,1);
    Enoise = sum(noise.^2);
    Esignal  = sum(yy.^2);
    
    K = Esignal * 10^(-SNR/10) ; % note factor of 10 because energies are squared
    noise_amp = sqrt(K / Enoise);
    
%     noise = noise * noise_amp;
%     10*log10(Esignal/sum(noise.^2))

    yy = yy + noise_amp*randn(size(yy));

end
