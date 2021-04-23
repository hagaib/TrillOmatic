function [harmonics , energies] = harmonics_of_signal(xx , fs , yin , maxharms)
%maxharms is the maximal number of harmonics desired (excluding f0)
%output: harmonics == matrix of detected harmonics in time instances
%specified in yin.time
f0=yin.f0;
dips = yin.dips;
time = yin.time;

winsize = 256;
harmonics=nan(length(f0) , maxharms);
energies=nan(length(f0) , maxharms+1);
for i=1:length(f0)
    if(dips(i) < 0.6)        
        indices = time(i)*fs;
        indices = max(1,floor(indices-winsize/2)):min(length(xx),floor(indices+winsize/2-1));
        [h , e] = find_harmonies(xx(indices) , fs , f0(i) , [-f0(i)/2 , f0(i)/2] , 2048);
        len = min(maxharms+1 , length(h));
        h(h==0) = nan;
        harmonics(i , 1:len-1)=h(2:len);
        energies(i , 1:len) = e(1:len);
    end
end
