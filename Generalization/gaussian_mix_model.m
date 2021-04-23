function model = gaussian_mix_model(x , points)
%GAUSSIAN_MIX_MODEL Summary of this function goes here
%   Detailed explanation goes here
if(mod(length(x),2) ~= 0)
    error('x must be of length of a multiple of 2');
end
ind = 1:2:length(x);
mu = x(ind);
sigma = x(ind+1);

model = zeros(length(points) , length(ind));
for i=1:length(ind)
    model(:,i) = exp(-0.5*((points - mu(i))/sigma(i)).^2);
end


