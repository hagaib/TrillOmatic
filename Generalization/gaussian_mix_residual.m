function residual = gaussian_mix_residual(x , points , known_val , isplot)
% residual function for nonlinear least squares curve fitting using
% 1D gaussian mixture model.
%inputs:
% x: model parameters [mu1, sigma1, mu2, sigma2, ... muN, sigmaN]
% points: range of data points (x values).
% known_val: value of data points (y values).
% is_plot: boolean. toggles plotting of optimization process for debugging
% purposes
%
%outputs:
% residual: residual to minimize (model - data)

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