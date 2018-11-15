
classdef DataExploration < handle
    
    properties
        
        stSubject;
        sFolderMain = pwd;
        sFolderFunctionsExternal = 'I:\IHAB_DataExtraction\functions';
        sFolderFunctions = 'functions';
        sFolderData = 'data';
        prefix;
        mColors;
        hSubPlot1;
        hSubPlot2;
        
    end
    
    methods
        
        function obj = DataExploration()
            
            addpath(genpath(obj.sFolderFunctions));
            addpath(genpath(obj.sFolderFunctionsExternal));
            
            if ismac
                obj.prefix = '/usr/local/bin/';
            else
                obj.prefix = '';
            end
            
            obj.stSubject = struct('Name', [], ...
                'Date', [], ...
                'Experimenter', [], ...
                'Folder', [], ...
                'Appendix', '-');
            
            obj.mColors = getColors();
            
            obj.analyse();
            
        end
        
        function [] = analyse(obj, ~, ~)
            
            obj.stSubject.Name = 'KE07IN22';
            obj.stSubject.Date = '180607';
            obj.stSubject.Experimenter = 'mh'; 
            obj.stSubject.Folder = [obj.sFolderMain, filesep, ...
                obj.sFolderData, filesep, obj.stSubject.Name, '_', ...
                obj.stSubject.Date, '_', obj.stSubject.Experimenter];
            obj.stSubject.Appendix = '';
            
    
            sFeatureFolder = [obj.stSubject.Folder, filesep, ...
                obj.stSubject.Name, '_AkuData'];
            
            
            tableDates = getdatesonesubject(obj);
            nDate = 1;
            desiredDay = tableDates.(obj.stSubject.Name)(nDate);
            desiredPart = 1;
            
            
            
            
%             tic
        szFeature = 'PSD';
        [DataPSD, timeVecPSD, ~] = getObjectiveDataOneDay(obj, desiredDay, szFeature, desiredPart);
%         toc
        %    save DataMat DataRMS timeVecRMS DataPSD timeVecPSD NrOfParts

%         tic
%         load('datapsd');
%         toc

        % Data conversion
        [Cxy,Pxx,Pyy] = get_psd(DataPSD);
%         clear DataPSD;
%         Cohe = Cxy./(sqrt(Pxx.*Pyy) + eps);


[nFrames, nFrameLength] = size(Pxx)
nStepSize = 128;
nLPCCoeffs = 12;
nFormants = 3;
nFs = 16000;
mFormants = zeros(nFrames, nFormants);
% Lower Limit in Hertz (only important for non-telephone voice)
nLowerBound = 200;
% Maximum Bandwidth
nBW = 400;

mSpec = zeros(nFrames, nFrameLength);

% hFig1 = figure();

for iFrame = 1:nFrames
     
    nIdxIn = (iFrame - 1) * nStepSize + 1;
    nIdxOut = nIdxIn + nFrameLength - 1;
    
    
    vSpec = Cxy(iFrame,:);
    
    vSpec(225:238) = mean([vSpec(204:224), vSpec(239:259)]);
    vSpec(335:355) = mean([vSpec(314:334), vSpec(356:376)]);
    vSpec(451:467) = mean([vSpec(430:450), vSpec(468:488)]);
    
    mSpec(iFrame, :) = vSpec;
    
%     if iFrame==1
%         hAx = plot(20*log10(abs(vSpec)));
%     else
%         hAx.YData = 20*log10(abs(vSpec));
%     end
%     drawnow;
    
    vSpecComplete = [vSpec(end:-1:1),vSpec(2:end-1)];
%     vSpecComplete = [vSpec(2:end-1),vSpec(end:-1:1)]
    
    vSample = ifft(circshift(real(vSpecComplete),nFrameLength));
   
    
%     vSample = vSignal(nIdxIn:nIdxOut);
    % LPC coefficients of sample
    vLPC = lpc(vSample,nLPCCoeffs);
    vLPC(isnan(vLPC)) = 0;
%     % Roots of the polynomial
    vRoots = roots(vLPC(1:3*nFormants));
%     % Only roots with positive imaginary part
    vRoots = vRoots(imag(vRoots)>=0);
%     % Sorted
    [vFormantsTemp,vIdx] = sort(angle(vRoots)/(2*pi)*nFs);
%     % Bandwidth
    vBW = -1/2*(nFs/(2*pi))*log(abs(vRoots(vIdx)));
%     % Ex-/Inclusion regarding features
    vFormantsTemp = trimData(vFormantsTemp,vBW,nLowerBound,nBW);
%     % Take the first three or less formants
    mFormants(iFrame,:) = guaranteeLength(vFormantsTemp,nFormants);
    
end
        save('formants','mFormants');

1;
        
        
        % OVD
%         stInfoOVDAlgo.fs = 16000;
%         [OVD_result_fixed, MeanCohere,~] = computeOVD_Coh(Cohe, timeVecPSD, stInfoOVDAlgo);
            
            
%             
%             
%             szFeature = 'RMS';
%             
%             
%             [DataRMS, timeVecRMS, ~] = getObjectiveDataOneDay(obj, ...
%                 desiredDay, szFeature, desiredPart);
%        
%             
% %             nInspectLength = 20000;
% %             DataRMS = DataRMS(1:nInspectLength,:);
% %             timeVecRMS = timeVecRMS(1:nInspectLength);
% 
%             hFig1 = figure();
%             obj.hSubPlot1 = subplot(2,1,1);
%             plot(timeVecRMS, 20*log10(abs(DataRMS)));
%             axis tight;
%             obj.hSubPlot1
%             
%             
%            
%             nFrameLength = 256;
%             nHopSize = 128;
%             nFrames = 1+ floor((length(DataRMS)-nFrameLength)/nHopSize);
% %             mFrames = zeros(nFrames, nFrameLength, 2);
%             mFramesCorr = zeros(nFrames, 2*nFrameLength-1);
%             vMaxLag = zeros(nFrames, 1);
%             vDateTimeCorr = linspace(timeVecRMS(1), timeVecRMS(end), nFrames);
% 
%             
%             for iFrame = 1:nFrames
%                
%                iIn = (iFrame - 1)*nHopSize + 1;
%                iOut = iIn + nFrameLength - 1;
%                
% %                mFrames(iFrame, :, :) = DataRMS(iIn: iOut, :);
%                [vCorr, vLag] = xcorr(DataRMS(iIn:iOut,1), DataRMS(iIn:iOut,2));
%                mFramesCorr(iFrame, :) = vCorr;
%                [nCorr, nIdx] = max(vCorr);
%                vMaxLag(iFrame) = log10(nCorr);%vLag(nIdx);
% 
%             end
%             
%             
%             obj.hSubPlot2 = subplot(2,1,2);
% %             hImage = imagesc(log10(mFramesCorr)');
% %             hImage.XData = datenum([vDateTimeCorr(1), vDateTimeCorr(end)]);
% %             hold on;
%             plot(vDateTimeCorr, vMaxLag);
% %             ylim([-nFrameLength,nFrameLength]);
% %             ylim([-1,1]);
%             
%             linkaxes([obj.hSubPlot1, obj.hSubPlot2], 'x');
           
            
        end
        
        function [] = moveAxis(obj, source, event)
           
            source
            event
            
        end
        
        
    end
    
end




% if ~isempty(timeVecRMS)
%     % No Inclusion of Parts shorter than e.g. 10 Minutes % UK
%     temp_time = timeVecRMS(end) - timeVecRMS(1);
%     temp_duration = duration(temp_time);
%     if (minutes(temp_duration) < obj.stPreferences.MinPartLength) && (obj.stPreferences.MinPartLength > 0)
%         return;
%     end
% end


%
% szFeature = 'PSD';
% [DataPSD, timeVecPSD, ~] = getObjectiveDataOneDay(obj, desiredDay, szFeature, desiredPart);
% %    save DataMat DataRMS timeVecRMS DataPSD timeVecPSD NrOfParts
%
%
%
%
% % Data conversion
% [Cxy,Pxx,Pyy] = get_psd(DataPSD);
% clear DataPSD;
% Cohe = Cxy./(sqrt(Pxx.*Pyy) + eps);
% if isempty(Cohe)
%     return;
% end
% % OVD
% stInfoOVDAlgo.fs = 16000;
% [OVD_result_fixed, MeanCohere,~] = computeOVD_Coh(Cohe, timeVecPSD, stInfoOVDAlgo);
%
% % OVD adaptive
% stInfoOVDAlgo.fs = 16000;
% stInfoOVDAlgo.adapThresh = 0.5;
% stInfoOVDAlgo.additive = p.Results.additive;
% OVD_result_adaptive = computeOVD_Coh(Cohe, timeVecPSD, stInfoOVDAlgo);

