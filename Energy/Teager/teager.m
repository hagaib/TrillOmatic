function yy = teager( xx )
% Computes Teager Energy Operator (TEO). For more information refer to  
% P. Maragos, J. F. Kaiser, and T. F. Quatieri, “Energy separation in 
% signal modulationswith application to speech analysis,”IEEE transactions 
% on signal processing, vol. 41,no. 10, pp. 3024–3051, 1993

   yy =  xx(2:end-1).^2 - xx(1:end-2).*xx(3:end);


