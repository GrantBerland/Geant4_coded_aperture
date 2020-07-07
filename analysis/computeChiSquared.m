function chiSquare = computeChiSquared(est, true, plotOn, figNum)

% This function should be used to compare the initial PDF and the
% fitted model that attempts to reconstruct that PDF, not the base
% reconstructed image.

%counting_variance = cov(est);
counting_variance = ones(256,256);

chiSquareMap = 1 ./ counting_variance .* (est - true).^2;

% Pearson chi-squared test
%chiSquareMap = 1./true .* (est - true).^2;

chiSquare = sum(sum( chiSquareMap ));

if plotOn == 1
    figure(figNum); subplot(1,2,2);
    contourf(chiSquareMap);
    colorbar(); title(sprintf("\\chi^2 Map, \\chi^2=%.3e", chiSquare));
end

end
