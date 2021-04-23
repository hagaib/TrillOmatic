function yy = prevent_audio_clipping(xx)
%PREVENT_CLIPPING Normalize signal xx so that all values lie in the
%interval [-1,1]. Essensial in preventing audio clipping when writing to
%wav format in 16 bit samples.
    
yy = xx;
a = min(xx);
b = max(xx); 
if(b>1 || a<-1)
    if(b-a > 2)
        yy = 2*(xx-a)/(b-a)-1;% copy entire range to [-1,1]
    elseif(b>1)
        yy = xx -b +1;%shift signal down by 1-b
    elseif(a<-1)
        yy = xx -a -1;%shift signal up by -1-a
    end
end
