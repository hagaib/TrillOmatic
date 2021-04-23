function yy = add_noise_to_signal(xx , snr , ignore0 , isCorrelated , noisetype)
% xx is assumed to be a clean signal
% snr is a vector of desired SNR values in db
% ignore0 - boolean optional parameter. If true , samples set to 0 are
% ignored in the calculation of mean energy in xx. Defaults to false
% isCorrelated - optional boolean. in case multiple SNRs are specified,
% determines if columns of output matrix yy are added with correlated white
% noise or uncorrelated white noise.
% noisetype - 'white' = AWGN , 'agmon' = synthetic agmon ha'hula noise
% output is new signal with additive white noise

% example: yy = add_noise_to_signal(xx , [5 , 0 , -5] , true , false , 'white');

if(nargin<3)
    ignore0 = false;
end

if(nargin<4)
    isCorrelated = false;
end

snr = snr(:)';


if(ignore0)
    xxe = mean(xx(xx~=0).^2);
else
    xxe = mean(xx.^2);
end

noisecoeffs = 10.^(-snr./20)*sqrt(xxe); 
% so that 10*log10(mean(yy.^2)/mean(noise.^2)) = snr

matsize = [length(xx) , length(snr)];
if(isCorrelated)
    matsize = size(xx);
end

if(isequal(noisetype , 'agmon'))
        noise = generate_agmon_noise(matsize(1) , matsize(2));
    else
        noise = randn(matsize); 
end
    
if(isCorrelated)
    %correlated noise:
    noise = repmat(noise , 1 , length(snr));
end

noisemat =  noise * diag(noisecoeffs);

yy = repmat(xx , 1 , length(snr));

yy = yy + noisemat;