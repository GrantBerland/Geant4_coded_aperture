function PSFest = computePSFestimate(est, true, RL_PSFest, plotOn)

% Minimum geometric resolution of coded aperture mask
%minPhysicalRes = 12.8; % pixels / 50 km
minPhysicalRes = 0; 

% x = [sigma_i, sigma_j] in model
% PSF_est = A * exp(-1/2 * (i^2/sigma_i^2 + j^2/sigma_j^2)
i = linspace(0, 255, 256); j = i;
[i,j] = meshgrid(i, j);

gaussian_kernel = @(x) 1/(2*pi*x(1)*x(2))*exp(-1/2*((i-128).^2/x(1)^2 + (j-128).^2/x(2)^2));

cost_fnc = @(x) sum(sum( ( conv2( true, gaussian_kernel(x) ,'same') - est ).^2 ));

% Initial guess at [sigma_x^2 , sigma_y^2] ~ U[min, 100]
x0 = rand(2,1)*(100-minPhysicalRes)+minPhysicalRes;

if plotOn == 1
    options = optimset('PlotFcns',@optimplotfval);
else
    options = optimset();
end

%PSFest = fminsearch(cost_fnc, x0, options);
%PSFest = fminunc(cost_fnc, x0, options);
PSFest = fmincon(cost_fnc, x0,[],[],[],[],...
                [minPhysicalRes;minPhysicalRes],[Inf;Inf],@mycon,options);

figure(); subplot(1,2,1);
contourf(gaussian_kernel(PSFest));
title("Gaussian Kernel PSF Estimate"); colorbar();

subplot(1,2,2);
contourf(RL_PSFest);
title("R-L Iterative PSF Estimate"); colorbar();
end
function [c,ceq] = mycon(x)
differenceTolerance = 10;
c   = abs(x(1) - x(2)) - differenceTolerance;
ceq = [];
end