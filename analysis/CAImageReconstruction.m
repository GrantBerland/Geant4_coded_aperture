classdef CAImageReconstruction
    % This class encapsulates the various coded aperture image
    % reconstruction methods as I implement them. The first come from the
    % work of Caroli et al. 1987 and are the direct balanced correlation
    % methods. Later methods include iterative maximum entropy methods
    % (MEMs) or minimum variance methods. 
    %
    % Written by Grant Berland, June 2020
    
    properties
        
        % Input members
        Filename;
        SNR;
        EnergyRange;
        SmoothingFactor;
        ReconstructionMethod;
        
        % Internal members
        RawImage;
        DarkImage;
        DeconvolvedImage;
        PSF_estimate;
        
        Decoder;
        Mask;
        
    end
    
    methods
        
        % Constructor
        function obj = CAImageReconstruction(inputArg1,inputArg2,inputArg3,inputArg4)
            
            if nargin == 4
                obj.Filename        = inputArg1;
                obj.SNR             = inputArg2;
                obj.EnergyRange     = inputArg3;
                obj.SmoothingFactor = inputArg4;
            else
                error("Please enter filename, SNR, energy range, and smoothing factor");
            end
            
            % Load mask and decoder
            tmpLoad     = load('./CA_files/decoder.mat');
            obj.Decoder = tmpLoad.decoder;
            tmpLoad     = load('./CA_files/decoder.mat');
            obj.Mask    = tmpLoad.mask;
            
            hitFileName    = sprintf("../data/hit_%s.csv", obj.Filename);    
            signalFileName = sprintf("../data/signal_%s.csv", obj.Filename);

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
            obj.RawImage   = zeros(16,16);
            obj.DarkImage  = zeros(16,16);
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
                    obj.RawImage = obj.RawImage + binnedData.Values + poissrnd(perDetRate/256, [16 16]);
                else
                    obj.DarkImage = binnedData.Values + poissrnd(perDetRate/256, [16 16]);
                end
            end
            close 1;

            % Normalize counts after image coaddition
            obj.RawImage = obj.RawImage / 11;
            
        end
        
        % Balanced correlation with G(null) := -1 
        % G(null) = tau/(tau-1), tau = M/N 
        % from [Caroli, et al. 1987]
        function outputArg = reconstructionMethod1(obj,inputArg)

            outputArg = obj.Property1 + inputArg;
            
            
        end
        
        % Balanced correlation with G(null) := -0.7059
        % G(null) = -lambda/(M - lambda), lambda ~ avg. sidelobe height
        % from [Caroli, et al. 1987]
        function outputArg = reconstructionMethod2(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
        
        
        
    end
    
end

