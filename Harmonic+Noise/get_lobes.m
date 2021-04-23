function [ peaks , ind_peaks , lobe_ampl, lobe_bounds] = get_lobes(P)
%find lobe information from Power Spetrum P!

[ext , ind] = extrema(P,0,1);
peaks = zeros(size(ext));
ind_peaks = zeros(size(ext));
lobe_ampl = zeros(size(ext));
lobe_bounds = zeros(length(ext),2);
k=1;
for i=2:(length(ind)-1)
    if(ext(i)>ext(i-1) && ext(i)>ext(i+1))
        %ext(i) is a peak
        peaks(k) = ext(i);
        ind_peaks(k) = ind(i);
        lobe = P(ind(i-1):(ind(i+1)-1));
        lobe_ampl(k) = sum(lobe);
        lobe_bounds(k,:) = [ind(i-1),(ind(i+1)-1)];
        k=k+1;
    end
end

peaks = peaks(1:k-1);
ind_peaks = ind_peaks(1:k-1);
lobe_ampl = lobe_ampl(1:k-1);
lobe_bounds = lobe_bounds(1:k-1,:);
