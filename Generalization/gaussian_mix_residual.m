function residual = gaussian_mix_residual(x , points , known_val , isplot)
%GAUSSIAN_MIX_MODEL Summary of this function goes here
%   Detailed explanation goes here

if(nargin<4) , isplot = false; end

model = gaussian_mix_model(x , points);


%% Find coefficients with linear least squares regression
global coeffs
coeffs = model \ known_val;
coeffs = max(coeffs , 1);
residual = model*coeffs - known_val;

if(isplot)
figure(isplot)
plot(points , known_val , points, model*coeffs)
title('Gaussian Model Fitting Optimization');
drawnow();
% coeffs
end