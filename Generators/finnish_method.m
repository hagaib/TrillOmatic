function [y , yindata] = finnish_method(x, fs , max_dur , yindata)
%FINNISH_METHOD Fixes track x with missing data using the finnish method.
% after "Reconstruction Method for Missing or Damaged Long Portions in 
% Audio Signal" by ISMO KAUPPINEN AND JYRKI KAUPPINEN
% max_length is the maximal interval length to fix (in seconds)

if(nargin < 3) , max_dur = 5 * 10^3; end % 2 ms
xorig = x;

s_extended_dur = 5;
s_synthesis_windur = floor(0.001 * fs);
s_synthesis_windur = s_synthesis_windur + 1 - mod(s_synthesis_windur , 2);

dx = x(2:end) - x(1:end-1);
indicator = medfilt1(double(dx==0) , 3);
segs_miss = logical2segments(indicator , fs);
segs_dur = segs_miss(2,:) - segs_miss(1,:);
b_segs = segs_dur < max_dur;
b_segs(1) = 0; b_segs(end) = 0;

% figure(1)
% plot(x) , hold on
% plot(indicator) , hold off
% legend('original' ,'reconstructed')

% x = sin(2*pi*440*(0:1/fs:8))';

% figure(2)
% plot(0,0);
% hold on

for j=1:length(b_segs)
if(~b_segs(j) ) , continue; end
s_start = floor( segs_miss(1 , j) * fs);
s_end = floor( segs_miss(2 , j) * fs);
seglen = s_end - s_start + 1 ;

s_leftseg = floor(fs*[segs_miss(2 ,j-1) , segs_miss(1 , j)]);
s_rightseg = floor(fs*[segs_miss(2 , j) , segs_miss(1 , j+1)]);
s_leftseg(2) = s_leftseg(2) - 1;
s_rightseg(1) = s_rightseg(1) + 1;
win = hanning((seglen+s_extended_dur)*2+1);
win = win(1:seglen+s_extended_dur);
synthwin = hanning(s_synthesis_windur);
synthwin = synthwin(1:floor(s_synthesis_windur/2));
segfixleft = zeros(seglen + s_extended_dur, 1);
segfixright = zeros(seglen + s_extended_dur , 1);
extrapolated_values = zeros(seglen + 2*s_extended_dur , 1);
lpcord = 20;%floor(seglen*0.05);

if(s_leftseg(2) - s_leftseg(1) < lpcord+s_extended_dur || ...
        s_rightseg(2) - s_rightseg(1) < lpcord+s_extended_dur)
    continue;
end

% plot(x) ,hold on;
x(s_leftseg(2)-length(synthwin)+1:s_leftseg(2)) = ...
    x(s_leftseg(2)-length(synthwin)+1:s_leftseg(2)) ./ wrev(synthwin);
x(s_rightseg(1):s_rightseg(1)+length(synthwin)-1) = ...
    x(s_rightseg(1):s_rightseg(1)+length(synthwin)-1) ./ synthwin;
% plot(x)

for i=1:seglen+s_extended_dur
    endseglpc = s_leftseg(2)+i-1-s_extended_dur;
    a = lpc(x(s_leftseg(1):endseglpc) , lpcord);
    segfixleft(i) = -a(2:end)*x(endseglpc:-1:endseglpc-lpcord+1);
    x(endseglpc+1) = segfixleft(i);
end
extrapolated_values(1:end-s_extended_dur) = segfixleft .* wrev(win);

% plot(s_start-s_extended_dur:s_end , segfixleft);


for i=1:seglen+s_extended_dur
    startseglpc = s_rightseg(1)-i+1+s_extended_dur;
    a = lpc(x(s_rightseg(2):-1:startseglpc) , lpcord);
    segfixright(i) = -a(2:end)*x(startseglpc:startseglpc+lpcord-1);
    x(startseglpc-1) = segfixright(i);
end
extrapolated_values(s_extended_dur+1:end) = ...
    extrapolated_values(s_extended_dur+1:end) + wrev(segfixright) .* win;

x(s_start-s_extended_dur:s_end+s_extended_dur) = extrapolated_values;


% plot(x), hold off

yind_start = find(yindata.time>=segs_miss(1,j) , 1);
yind_end = find(yindata.time<=segs_miss(2,j) , 1 ,'last');
yinan = isnan(yindata.f0(yind_start:yind_end));
while(sum(isnan(yindata.f0(yind_start:yind_end))))
    indyinan = yind_start + find(yinan==1,1) - 1;
    left_indyinan = indyinan;
    right_indyinan = indyinan;
    while(isnan(yindata.f0(left_indyinan))) , left_indyinan = left_indyinan-1; end
    while(isnan(yindata.f0(right_indyinan))) , right_indyinan = right_indyinan+1; end
    yindata.f0(left_indyinan:right_indyinan) = ...
        linspace(yindata.f0(left_indyinan) , yindata.f0(right_indyinan) , right_indyinan-left_indyinan+1);
    yinan = isnan(yindata.f0(yind_start:yind_end));
end

end

y = x;


end