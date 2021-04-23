function yy = teager( xx )
%TEAGER Summary of this function goes here
%   Detailed explanation goes here

   yy =  xx(2:end-1).^2 - xx(1:end-2).*xx(3:end);


