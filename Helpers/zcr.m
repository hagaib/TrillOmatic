function y=zcr(xn,Ts,flag)
% Zero Crossing rate
% flag = 0: counts zero crossings
% flag = 1: calculates zero crossing rate (in herz)
if flag==0
y=sum(sign(xn(1:end-1).*xn(2:end))<1);
else
    dur=length(xn)*Ts;
    y=sum(sign(xn(1:end-1).*xn(2:end))<1)/(2*dur);
end
    
