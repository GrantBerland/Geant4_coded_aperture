
%{
% To run on work computer:
%
% Send photon and run files to work computer
%
% FILENAME='test3';
%
% scp ./photonFiles/$FILENAME_photons.csv ./build/runFile_$FILENAME.mac gdberla@lair13.int.colorado.edu:/home/gdberla/Documents/Research/AEPEX/shielding_simulation/analysis/photonFiles
% scp ./build/runFile_$FILENAME.mac gdberla@lair13.int.colorado.edu:/home/gdberla/Documents/Research/AEPEX/shielding_simulation/build
% 
% Run simulation
%
% ssh gdberla@lair13.int.colorado.edu && cd /home/gdberla/Documents/Research/AEPEX/shielding_simulation/build && ./main runFile_$FILENAME.mac > output.log
% 
% Return results
%
% scp gdberla@lair13.int.colorado.edu:/home/gdberla/Documents/Research/AEPEX/shielding_simulation/data/{signal,hit}_$FILENAME.csv ../data/
%}

close all;

%%%%%% User inputs %%%%%%

fileName          = 'test3';  % No file extension
number_of_photons = 100;   % [photons]
photonsPerSource  = 100;
photon_powerLaw_k = 2.704; % []
SNR               = 25;   % []

smoothingFactor    = 2;  % 1 = no smoothing
chiSquareMapPlotOn = 1;
optimizationPlotOn = 0;

%k(100 keV) = 3.258
%k(200 keV) = 2.874
%k(300 keV) = 2.704

%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a drawing, estimate a PDF, and sample point sources to simulate
generateImage(fileName, number_of_photons, photonsPerSource);

% Write simulation macro file and run simulation using the generated photons
runSimulation(number_of_photons*photonsPerSource, photon_powerLaw_k, fileName);

% Image reconstruction algorithm
[image_estimate, RL_PSFest] = reconstructImageMethod1(fileName, SNR, [50 300], smoothingFactor);

% Plot original estimate in convolution space
figure(1);
contourf(image_estimate);
colorbar();
title("Original Image Estimate");

% Resize image to cast on to 256x256 image 
image_estimate = resizem(image_estimate, 256/63);

% Retrieve original drawn image PDF for comparison
[z, fig] = getOriginalImageData(fileName);

figure(fig); subplot(1,2,2); hold on;
contour(image_estimate, 'DisplayName','Reconstructed Image');
title("Image Estimate");
colorbar();

% Take the difference of the normalized reconstruction and PDF
residuals = image_estimate / sum(sum(image_estimate)) - z; 

errorFigNum = 3;
chiSquare = computeChiSquared(image_estimate, z, chiSquareMapPlotOn, errorFigNum);
fprintf("Chi squared score: %.3e\n", chiSquare);

figure(errorFigNum);
subplot(1,2,1);
contourf(residuals);
title("Reconstruction Residuals");
colorbar();

%PSFest = computePSFestimate(image_estimate, z, ...
%                resizem(RL_PSFest, 256/33), optimizationPlotOn);
%fprintf("sigma_x^2 = %.3e , sigma_y^2 = %.3e\n", [PSFest(1), PSFest(2)]);

% Optional: visualize generated photons
%ph = downsample(importfile_photonFile(fileName), photonsPerSource);
%quiver3(ph.x, ph.y, ph.z, ph.px, ph.py, ph.pz, 20)

function generateImage(fileName, Nsamples, photonsPerSource)
% A distribution or shape is chosen to be generated in an array of 
% arbitrary resolution

function [xq,yq,z] = computeGrid(x1,x2,fout,res)
[xq,yq] = meshgrid(linspace(min(x1),max(x1), res) ,...
                   linspace(min(x2),max(x2), res) );
orig_state = warning;
warning('off','all');
z = griddata(x1,x2,fout,xq,yq);
warning(orig_state);
end

imtool close all;  % Close all figure windows created by imtool.
workspace;  % Make sure the workspace panel is showing.
figure(1); clf; subplot(1,2,1);
fontSize = 12;

grayImage = zeros([256,256]);
imshow(grayImage, []); hold on;
t = linspace(0, 2*pi);
xCircle = 128 * (cos(t) + 1);
yCircle = 128 * (sin(t) + 1);
plot(xCircle, yCircle); hold off;

set(gca, 'ydir', 'normal');
axis on;
title('Draw Here', 'FontSize', fontSize);
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
message = sprintf("Left click and hold to begin drawing a freehand path.\nSimply lift the mouse button to finish.");
uiwait(msgbox(message));
% User draws curve on image here.
hFH = imfreehand();
% Get the xy coordinates of where they drew.
xy = hFH.getPosition;
% get rid of imfreehand remnant.
delete(hFH);
% Overlay what they drew onto the image.
hold on; % Keep image, and direction of y axis.
xCoordinates = xy(:, 1);
yCoordinates = xy(:, 2);
plot(xCoordinates, yCoordinates, 'ro', 'LineWidth', 2, 'MarkerSize', 5);
% Ask user if they want to burn the line into the image.
promptMessage = sprintf('Is this the shape you want to generate?');
titleBarCaption = 'Continue?';
button = questdlg(promptMessage, titleBarCaption, 'Yes', 'No', 'Yes');
if strcmpi(button, 'Yes')
  cla;
  hold off;
  for k = 1 : length(xCoordinates)
    row = int32(yCoordinates(k));
    column = int32(xCoordinates(k));
    row(row <= 0) = 1; row(row > 256) = 256;
    column(column <= 0) = 1; column(column > 256) = 256;
    grayImage(row, column) = 255;
  end
  imshow(grayImage, []); set(gca, 'ydir', 'normal');
  axis on;
  caption = sprintf("EPP Drawing (it's wonderful)");
  title(caption, 'FontSize', fontSize);

elseif strcmpi(button, 'No')
    % Retry
    generateImage(fileName);
end

xCoordinates = xy(:, 1);
yCoordinates = xy(:, 2);
numberOfKnots = length(xCoordinates);
% Close gaps that you get when you draw too fast.
% Use splines to interpolate a smoother curve,
% with 10 times as many points,
% that goes exactly through the same data points.
samplingRateIncrease = 10;
newXSamplePoints = linspace(1, numberOfKnots, numberOfKnots * samplingRateIncrease);
% smoothedY = spline(xCoordinates, yCoordinates, newXSamplePoints);
% Make the 2D array where the top row is the x coordinates and the bottom row is the y coordinates,
% but with the exception that the left column and right column is a vector that gives the direction of the slope.
yy = [0, xCoordinates', 0; 1, yCoordinates', 1];
pp = spline(1:numberOfKnots, yy); % Get interpolant
smoothedY = ppval(pp, newXSamplePoints); % Get smoothed y values in the "gaps".
% smoothedY is a 2D array with the x coordinates in the top row and the y coordinates in the bottom row.
smoothedXCoordinates = smoothedY(1, :);
smoothedYCoordinates = smoothedY(2, :);
% Plot smoothedY and show how the line is
% smooth, and has no sharp bends.
hold on; % Don't destroy the first curve we plotted.
hGreenCurve = plot(smoothedXCoordinates, smoothedYCoordinates, '-g');
title("EPP Drawing (it's wonderful)"', 'FontSize', fontSize);
% But smoothedXCoordinates and smoothedYCoordinates are not in pixel coordinates, they have fractional values.
% If you want integer pixel values, you have to round.
intSmoothedXCoordinates = int32(smoothedXCoordinates);
intSmoothedYCoordinates = int32(smoothedYCoordinates);
% But now it's possible that some coordinates will be on the same pixel if that's
% how they rounded according to how they were located to the nearest integer pixel location.
% So use diff() to remove elements that have the same x and y values.
diffX = [1, diff(intSmoothedXCoordinates)];
diffY = [1, diff(intSmoothedYCoordinates)];
% Find out where both have zero difference from the prior point.
bothZero = (diffX==0) & (diffY == 0);
% Remove those from the arrays.
finalX = intSmoothedXCoordinates(~bothZero);
finalY = intSmoothedYCoordinates(~bothZero);
% Now remove the green line.
delete(hGreenCurve);
% Plot the final coordinates.
plot(finalX, finalY, '-y');

resolution = 500;
gridx1 = linspace(0, size(grayImage, 1), resolution);
gridx2 = linspace(0, size(grayImage, 2), resolution);
[x1,x2] = meshgrid(gridx1, gridx2);
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];

subplot(1,2,2); hold on; grid on; axis equal;
title('EPP Source PDF', 'FontSize', fontSize);

% Returns a kernel density estimate of the drawing
den = ksdensity([finalX; finalY]',xi, 'Bandwidth', 10);

[xq,yq,z] = computeGrid(xi(:,1),xi(:,2),den,resolution);
contourf(xq,yq,z);
xlabel('x_1'); ylabel('x_2');

% Create 2D CDF from PDF 
%zCDF = cumsum(cumsum(z, 1), 2);

% Rejection sampling to obtain N random samples
samples = zeros(Nsamples, 2);
Ntrial_samples = Nsamples*1000; % Memory versus speed tradeoff here
current_ind = 1;
while current_ind < Nsamples
    
    % Trial values to be tested, (U[1,256], U[1,256])
    vals = ceil(rand(Ntrial_samples, 2)*(size(z,1)-1) + 1);
    
    % Uniform randoms to compare with trial sample points
    comparison_rands = rand(Ntrial_samples,1);
    
    % Convert trial samples to linear indices
    ind_in = sub2ind(size(z), vals(:,1), vals(:,2));
    
    % Compare and accept if U[0,1] < f_XY(trial_X, trial_Y)
    accepted_idx = comparison_rands < z(ind_in);
    
    % Append accepted samples to array
    tmp = vals(accepted_idx,:);
    if length(tmp) > Nsamples
        tmp = tmp(1:Nsamples,:);
    end
    samples(current_ind:current_ind+sum(accepted_idx)-1,:) = tmp;
    if sum(accepted_idx) ~= 0
        current_ind = current_ind + sum(accepted_idx)-1;
    end
end

% Take only the first Nsamples accepted samples
if length(samples) > Nsamples
    samples = samples(1:Nsamples, :);
end

% Plot accepted samples
scatter(xq(1,samples(:,2)),yq(samples(:,1),1), 50, 'r.');
legend('f_{X_1,X_2}(x_1,x_2)','Samples');

boxSize = 30; % cm x cm

convFactor = 45/128; % deg/pixel
shiftFactor = 127;

% Shift from [0,255] px. -> [-127,128] px. -> (-45,45] deg
xSamples = (xq(1,samples(:,1))' - shiftFactor) * convFactor;
zSamples = (yq(samples(:,2),1)  - shiftFactor) * convFactor;


theta = atan2d(zSamples, xSamples);
phi   = sqrt(xSamples.^2 + zSamples.^2);


% Spherical definition of (x,y,z) from (r = unity, theta, phi)
x_comp =  sind(phi).*cosd(theta);
z_comp =  sind(phi).*sind(theta);
y_comp = -cosd(phi);

% Reproduce the Nsamples points photonsPerSource number of times
x_comp = reshape(repmat(x_comp, 1, photonsPerSource), length(x_comp)*photonsPerSource, 1);
z_comp = reshape(repmat(z_comp, 1, photonsPerSource), length(z_comp)*photonsPerSource, 1);
y_comp = reshape(repmat(y_comp, 1, photonsPerSource), length(y_comp)*photonsPerSource, 1);

% Uniformly random samples on the square domain, 
% constant y height above detector
x = (rand(Nsamples*photonsPerSource, 1)*2-1)*boxSize/2;
z = (rand(Nsamples*photonsPerSource, 1)*2-1)*boxSize/2;
y = 10*ones(Nsamples*photonsPerSource,1);

% Saves drawing and PDF figure
saveString = ['./originalImages/', fileName, '.fig'];
saveas(1, saveString);

% Write sampled photons to file
csvwrite(['./photonFiles/',fileName,'_photons.csv'], [x, z, -y, x_comp, z_comp, -y_comp]);
 
end
function runSimulation(number_of_photons, photon_E0, fileName)

runFileName = ['runFile_', fileName, '.mac'];

runFileString = {"# Set low energy EM processes, including fluorescence"
"/process/em/fluo true"
"/process/em/auger false"
"/process/em/augerCascade false"
"/process/em/pixe true"
"/process/em/deexcitationIgnoreCut false"
""
"/run/initialize"
""
"/control/verbose 0"
"/run/verbose 0"
"/event/verbose 0"
"/tracking/verbose 0"
""
"# Set hit and signal file names"
sprintf("/dataCollection/setHitFileName signal_%s.csv", fileName) 
""
sprintf("/energy/setPhotonFileName %s_photons.csv", fileName)
""
"# E0 of background distribution out of {100,200,300} keV"
sprintf("/energy/setFoldingEnergy %.3f # [keV]", photon_E0) 
""
"/energy/setDistributionType 5"
""
sprintf("/run/beamOn %d", number_of_photons)
""};

% Write to file
fID = fopen(['../build/', runFileName], "w");
fprintf(fID, '%s\n', runFileString{:});
fclose(fID);

fprintf("Running simulation with %s using %d photons...", runFileName, number_of_photons)

% Switch to simulation directory
cd ../build/

% Run simulation with runFileName
system(sprintf("./main %s > output.log", runFileName));

% Switch back to this directory
cd ../analysis/

fprintf("simulation complete!\n")

end
function [image_est, PSF_est] = reconstructImageMethod1(fileName, SNR, energyRange, smoothingFactor)

lowE  = energyRange(1);
highE = energyRange(2);

darkDetectorID = 14;

hitFileName    = sprintf("../data/hit_%s.csv", fileName);    
signalFileName = sprintf("../data/signal_%s.csv", fileName);

data = importfile_resultsFile(signalFileName);
hits = importfile_resultsFile(hitFileName);

detectorIDs = unique([data.det; hits.det]);

Nsig       = size(data, 1); % c/s/intrument, total signal photons 
Nintrinsic = 20; % c/s/det, intrinsic detector noise
dt         = 10; % sec, integration time
ndet       = 11; % number of detectors

perDetRate = Nsig^2 * dt / SNR^2 - Nsig - Nintrinsic * ndet;

if perDetRate < 0; perDetRate = 0; end

perDetRate = 10;

pixelArray     = cell(length(detectorIDs), 2);
rawIm          = zeros(16,16);
darkIm         = zeros(16,16);
figure(1);
for detector = 1:length(detectorIDs)
    
    iPixels = [data.i(data.det == detectorIDs(detector) & data.E > lowE & data.E < highE) ; ...
        hits.i(hits.det == detectorIDs(detector) & hits.E > lowE & hits.E < highE)];
    
    jPixels = [data.j(data.det == detectorIDs(detector) & data.E > lowE & data.E < highE) ; ...
        hits.j(hits.det == detectorIDs(detector) & hits.E > lowE & hits.E < highE)];
   
    binnedData = histogram2(iPixels, jPixels, 0:1:16, 0:1:16);
    
    pixelArray{detector,1} = detectorIDs(detector);
    pixelArray{detector,2} = binnedData.Values;
    
    if detectorIDs(detector) ~= darkDetectorID 
        rawIm = rawIm + binnedData.Values + poissrnd(perDetRate/256, [16 16]);
    else
        darkIm = binnedData.Values + poissrnd(perDetRate/256, [16 16]);
    end
end
close 1;

% Normalize counts after image coaddition
rawIm = rawIm / 4;

load('CA_files/decoder.mat','decoder');
load('CA_files/mask.mat'   ,'mask');

% This value ensures that the sum over the image is equal to 0
% G = -1

% This value ensures that the sum over the image is equal to the counts
% over the detector. Lambda is the average height of the sidelobes of the
% autocorrelation of the mask with itself.
%
% G = -lambda / (M - lambda) = -60/(145-60)
decoder(decoder == -1) = -0.7059;  

decoder = rot90(decoder, -1);

% Test value
%decoder(decoder == -1) = -0.4;

avgSigCounts = sum(sum(rawIm));
backgroundCounts  = sum(sum(darkIm));

totalFlux = 4.656 * (avgSigCounts - backgroundCounts);

naivePSF = conv2(mask, decoder);

rawIm = rawIm - darkIm;

%rawIm = rawIm - mean(rawIm);
rawDeconv = conv2(rawIm, decoder);
%rawDeconv(rawDeconv < 0) = 0;

rawDeconv = interp2(rawDeconv, 1.1);
%[image_est, PSF_est] = deconvblind(rawDeconv, naivePSF);

% Bypass iterative Lucy-Richardson deconvolution
image_est = rawDeconv;
PSF_est   = naivePSF;


% Normalize s.t. sum across image is 1
image_est = image_est / sum(sum(image_est));

% Assign signal flux to image and correct for inversion
image_est = rot90(image_est * totalFlux, 3);

image_est = imgaussfilt(image_est, smoothingFactor);

end
function [z, fig] = getOriginalImageData(fileName)

% Open original image and get data back from it
fig = openfig(['./originalImages/' fileName '.fig']);

dataObjs = findobj(fig,'-property','ZData');
z = resizem(dataObjs(2).ZData, 256/500);

end
