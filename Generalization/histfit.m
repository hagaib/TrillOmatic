
fun = @(x) gaussian_mix_residual(x , F, fhist);
x0 =   [1000 1000, 2000 1000]';
options = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt');
[x,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,[],[],options);

% global coeffs;
% figure(10)
% plot(F , fhist) , hold on;
% plot(F , gaussian_mix_model(x , F)*coeffs); hold off

figure(12)
f = fit(F, fhist ,'gauss2');
plot(f,F,fhist)

