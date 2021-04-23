function [harmonies , energies] = find_harmonies(xx , fs , f0 , search_interval , fft_size)
% finds maximal harmony based on peak analysis of fft of windowed signal xx
% length of xx should be even

ampdiff_threshold = 8; %db
spectenergy_ratio_threshold = 2;
freq_relative_error_threshold = 0.1;

len = fft_size;
bin_size = fs/len;
f = abs(fft(xx.*hamming(length(xx)) , fft_size));

if (mod(length(f),2))
    P = f(1:ceil(len/2));
    P(2:end) = P(2:end) + flipud(f(ceil(len/2)+1:end));
else
    P = f(1:len/2+1);
    P(2:end-1) = P(2:end-1) + flipud(f(len/2+2:end));
end

total_spectral_energy = sum(P);

[peaks , ind_peaks , lobe_ampl] = get_lobes(P);

% [ext , ind] = extrema(P,0,1);
% peaks = zeros(size(ext));
% ind_peaks = zeros(size(ext));
% lobe_ampl = zeros(size(ext));
% k=1;
% for i=2:(length(ind)-1)
%     if(ext(i)>ext(i-1) && ext(i)>ext(i+1))
%         %ext(i) is a peak
%         peaks(k) = ext(i);
%         ind_peaks(k) = ind(i);
%         lobe = P(ind(i-1):(ind(i+1)-1));
%         lobe_ampl(k) = sum(lobe);
%         k=k+1;
%     end
% end
% 
% peaks = peaks(1:k-1);
% ind_peaks = ind_peaks(1:k-1);
% lobe_ampl = lobe_ampl(1:k-1);

l=2;
harmonies = ones(size(peaks))*(-fs);
harmonies(1) = f0;
energies = zeros(size(harmonies));
energies(1) = floor(f0/bin_size) +1;

while(abs(harmonies(l-1))+f0 < fs/2)
    range = abs(harmonies(l-1))+f0 + search_interval;
    range_bins = floor(range/bin_size) +1;
    range_bins(2) = min(length(P) , range_bins(2));
    
%     plot(P(range_bins(1):range_bins(2)))
    
    peaks_in_range = find(ind_peaks>range_bins(1) & ind_peaks<range_bins(2));
    
    if (isempty(peaks_in_range))
        harmonies(l) = -l*f0;
        l=l+1;
        continue
    elseif(length(peaks_in_range)==1)
        spect_energy = sum(P(range_bins(1):range_bins(2)));
        lobe_energy = lobe_ampl(peaks_in_range);
        if(lobe_energy/(spect_energy-lobe_energy)>spectenergy_ratio_threshold && ...
                lobe_energy/total_spectral_energy > 0.1)
            harmonies(l) = (ind_peaks(peaks_in_range)-1)*bin_size;
            energies(l) = lobe_energy;
        else
            harmonies(l) = -l*f0;
        end
        l=l+1;
        continue
    end
    
    [w , ind_w] = max(peaks(peaks_in_range));
    tmp_range = peaks_in_range;
    tmp_range(ind_w)=[];
    w2 = max(peaks(tmp_range));
    amp_diff = 20*log10(w/w2);
    
    ind_w = peaks_in_range(ind_w);
    mean_spect_energy_in_range = (sum(lobe_ampl(peaks_in_range)) - lobe_ampl(ind_w))/(length(peaks_in_range)-1);
    
    h = (ind_peaks(ind_w)-1)*bin_size;
    if(     lobe_ampl(ind_w)/total_spectral_energy > 0.015 && ...%             amp_diff>ampdiff_threshold && ...
            lobe_ampl(ind_w)/mean_spect_energy_in_range > spectenergy_ratio_threshold && ...
            abs(h - l*f0)/(l*f0) < freq_relative_error_threshold && ...
            abs(h - l*f0) < 500  )
        harmonies(l) = h;
        energies(l) = lobe_ampl(ind_w);
    else
        if(abs(h - l*f0)/(l*f0) > freq_relative_error_threshold)
            h = l*f0;
        end
        harmonies(l) = -h;
    end

    l=l+1;
end

while (harmonies(l-1) < 0 || isnan(harmonies(l-1)))
    l=l-1;
end

harmonies = harmonies(1:l-1);
harmonies(harmonies<0) = 0;
energies = energies(1:l-1);
energies = 20*log10(energies);