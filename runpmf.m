%% Define the PMF type

pmfType = 'HR'; % type of the PMF ('UMR', 'bin', or 'HR')

%% Define the data source

% time resolution (aka pre-averaging time) of the data
aveTimeMins=30/60; % (min)

% if using data, mz, and t variables
dataSource=struct(...
    'mzs_bin',[],...
    'mzs_UMR',[],...
    'mzs_HR',mz,... % mzs_ parameters are 1xn-sized vectors (in Th)
    'dataMatrix_bin',[],...
    'dataMatrix_UMR',[],...
    'dataMatrix_HR',data,... %dataMatrix_ parameters are mxn-sized matrices (in cps)
    'compoundNames',{compoundNames},... % compoundNames is a 1xn-sized cell array containing the names of the HR-ions
    't',t); % t is a 1xm-sized datenum vector containing the time stamps

% if using a file preprocessed in tofTools (note that tofTools package is needed to load the file)
% dataSource='./data/paramEval1.h5';

%% Define the time range for the data to be analyzed

% kx2-sized datetime matrix listing the start and end times of the k good time periods
goodTimes=[datetime(2000,1,1,0,0,0) datetime(2100,1,1,0,0,0);];

%% Define global parameters for PMF 

initializeVariables = 1; % Running this code the first time, put 1 for this 
% so that the required variables are initialized (takes some time). 
% If global parameters or pmfType are modified, initialization has to be re-run.

useGPU = 0; % if GPU is used to boost the computation

mzRange = [40:45 47:54 58:72 74:400]; % what masses are included in the UMR-PMF or binPMF runs (always listed as UMR-type masses)
binWidth_bin = 0.01; % width of a bin in binPMF (in Th)
binWidth_UMR = nan; % width of a "bin" in UMR-PMF (in Th); should normally be 1, if less it's handled as bins, if nan it's calculated from lower and upper mass defects
lowerMassDefect = 0.25; % UMR minus this value (in Th) is the lower limit of a mass
upperMassDefect = 0.35; % UMR plus this value (in Th) is the upper limit of a mass

hrMassesToSkip = {'NO+','NO2+'}; % masses that will be skipped in HR-PMF run

% errorParameters for UMR, bin, and HR PMF separately
% The code assumes a following function to calculate the error matrix (in cps):
% sqrt(errorParameter1^2*dataMatrix+errorParameter2^2)
% For the first PMF run, the default values can be used for a quick PMF
% run. Then, the necessary variables are loaded for a given pmfType, and
% you can run errorParameterCalculation.m to derive the case-specific and
% pmfType-specific parameters.
% Note that the parameters depend on the instrument, time-averaging, etc.
errorParameter1_UMR=0.52;
errorParameter2_UMR=0.3;
errorParameter1_bin=0.52;
errorParameter2_bin=0.3;
errorParameter1_HR=1.31;
errorParameter2_HR=0.66;

%% Define parameters for the PMF run

% If any factor(s) is/are to be constrained, define the number of the
% constrained factors and their time series (TS) and mass spectra (MS).
% Note that both TS and MS have to be defined.
nConstrainedFactors = 0; % set to 0 for non-constrained run
constrainedTS = nan; % if constrained, define a mxp-sized matrix, where p is the number of constrained factors
constrainedMS = nan; % if constrained, define a pxn-sized matrix, where p is the number of constrained factors
involveConstrainedFactorsInQCalc = 1; % set to 0 for not including the constrained factors in the Q calculation in a constrained PMF run

% How many times PMF is repeated using a random seed. Initially it could be
% 10-50, but later it could be 100-500 after the final number of factors is decided.
repetitions = 20; 

% The range of number of factors. 
% Initially it could be 2:15, but later it should be fixed to a single value (after deciding the final number of factors).
pRange = 2:6; 

% Type of the PMF computation algorithm ('ALS', 'MULT', 'ALS5MULT1', or 'PSPMF')
pmfAlg = 'ALS5MULT1'; % default is 'ALS5MULT1', which means running alternating least squares (ALS) 5 times and then multiplicative solver (MULT) once

% PSPMF-specific settings
pspmfRepeats = {[5 1], [5 1], 1}; % numbers of repeats for algorithm inside PSPMF: {UMR normal time resolution run, bin lower time resolution run, bin normal time resolution run} 
pspmfAlgs = {'AM', 'AM', 'M'}; % algorithms inside PSPMF: {UMR normal time resolution run, bin lower time resolution run, bin normal time resolution run} (A means ALS, M means MULT)
timeResolutionLoweringFactor = 20; % How much the time are time-averaged for PSPMF. The default 20 means that the data are averaged into 20-point chunks.

% Plotting parameters
plotUnexpAsTime = 1; % 1 if unexplained variation is plotted as time series, 0 if as mass spectra
addMassTextDataTips = 1; % if labels of the peaks are added as data tips in HR-PMF plotting
printMassTexts = 1; % if labels of the tallest peaks are shown in HR-PMF plotting
printOnlyHighestMassTexts = 3; % the number of the tallest peaks' labels shown in HR-PMF plotting

% special parameters
downWeightingBinOn = 1; % if downweighting is enabled for binPMF
robustModeOn = 0; % set to 1, if abs 4 robust parameter is used also in the final calculation (internally it's always on)
umrQQexpTrackedWithBinPMF = 0; % advanced parameter for PSPMF to enable if Q is tracked for UMR instead or bin data

%% Running PMF

if initializeVariables %initializing required variables
    hw=waitbar(0,'Loading data');
    close(hw)

    if ischar(dataSource) % initializing data from tofTools h5 file

        dataMatrix_UMR=[];
        dataMatrix_bin=[];
        dataMatrix_HR=[];

        if ~strcmp(pmfType,'HR') % UMR or bin
                    
            mzs_bin=[];
            mzs_UMR=[];
                
            lowerBins=floor(lowerMassDefect/binWidth_bin);
            upperBins=floor(upperMassDefect/binWidth_bin);
                
            for i=1:length(mzRange)
                mzs_bin=[mzs_bin mzRange(i)];
        
                for j=1:lowerBins
                    mzs_bin=[mzs_bin mzRange(i)-j*binWidth_bin];
                end
        
                for j=1:upperBins-1
                    mzs_bin=[mzs_bin mzRange(i)+j*binWidth_bin];
                end
        
            end
        
            mzs_bin=sort(mzs_bin);
        
            if binWidth_UMR==1
                 mzs_UMR=mzRange-binWidth_UMR/2; %the normal calculation (using full range)
            elseif isnan(binWidth_UMR) % using only partial range
                 mzs_UMR=mzRange-lowerMassDefect;
                 binWidth_UMR=upperMassDefect+lowerMassDefect;
            else

                %the special calculation below
                lowerBins_UMR=floor(lowerMassDefect/binWidth_UMR);
                upperBins_UMR=floor(upperMassDefect/binWidth_UMR);
        
                for i=1:length(mzRange)
                    mzs_UMR=[mzs_UMR mzRange(i)];
        
                    for j=1:lowerBins_UMR
                        mzs_UMR=[mzs_UMR mzRange(i)-j*binWidth_UMR];
                    end
        
                    for j=1:upperBins_UMR-1
                        mzs_UMR=[mzs_UMR mzRange(i)+j*binWidth_UMR];
                    end
        
                end
        
                mzs_UMR=sort(mzs_UMR);
            end
            
        
            if exist('tof_get_convcoef.m','file')==0
                error('tofTools package is needed to initialize variables from a tofTools h5 file')
            end

            convcoef=tof_get_convcoef(dataSource);
            [pars,eq] = tof_get_mz_calib_pars(dataSource);
            funcArgs = eq();
                    
            clc
        
            disp('Importing UMR and bin data')
            
            if useGPU
                dataMatrix_UMR=gpuArray(single([]));
                dataMatrix_bin=gpuArray(single([]));
            else
                dataMatrix_UMR=single([]);
                dataMatrix_bin=single([]);
            end
                
           
            linestoread=1000;
        
            datasplitting1=1:linestoread:h5info(dataSource,'/tofData/signal').Dataspace.Size(2);
            datasplitting2=datasplitting1+linestoread-1;
            datasplitting2(end)=h5info(dataSource,'/tofData/signal').Dataspace.Size(2);
            toflength=h5info(dataSource,'/tofData/signal').Dataspace.Size(1);
            
            for iPart=1:length(datasplitting1)
                hw2=waitbar(iPart/length(datasplitting1),'Importing UMR and bin data');
                
                if iPart==length(datasplitting1)
                    linestoread=datasplitting2(end)-datasplitting1(end)+1;
                end
                
                tofsignal=h5read(dataSource,'/tofData/signal',[1 datasplitting1(iPart)],[toflength linestoread]); 
                tofsignalPart_UMR=zeros(size(tofsignal,2),length(mzs_UMR));
                tofsignalPart=zeros(size(tofsignal,2),length(mzs_bin));
        
                convcoefPart=convcoef(datasplitting1(iPart):datasplitting2(iPart));
                parsPart=pars(:,datasplitting1(iPart):datasplitting2(iPart));
        
                for i=1:size(tofsignalPart_UMR,2)
                    for j=1:size(tofsignalPart_UMR,1)
                        startI=funcArgs.inverseEQ(mzs_UMR(i),parsPart(:,j));
                        endI=funcArgs.inverseEQ(mzs_UMR(i)+binWidth_UMR,parsPart(:,j));
        
                        if startI>0&&endI>0
                            tofsignalPart_UMR(j,i)=sum(tofsignal(round(startI):round(endI),j))*convcoefPart(j);
                        end
                    end
                end
                
                for i=1:size(tofsignalPart,2) %masses
                    for j=1:size(tofsignalPart,1) %times
                        startI=funcArgs.inverseEQ(mzs_bin(i),parsPart(:,j));
                        endI=funcArgs.inverseEQ(mzs_bin(i)+binWidth_bin,parsPart(:,j));
                        
                        startI=startI+0.5;
                        endI=endI+0.5;
        
                        if startI>0&&endI>0
                            tofsignalPart(j,i)=(sum(tofsignal(floor(startI):(ceil(endI)-1),j))+...
                                (floor(startI)-startI)*tofsignal(floor(startI),j)+...
                                (endI-ceil(endI))*tofsignal(ceil(endI)-1,j))*convcoefPart(j);
                        end
                    end
                end
        
                dataMatrix_UMR=[dataMatrix_UMR; tofsignalPart_UMR];
                dataMatrix_bin=[dataMatrix_bin; tofsignalPart];
        
            end
        
            if useGPU
                mzs_bin=gpuArray(mzs_bin);
                mzs_UMR=gpuArray(mzs_UMR);
            end
        
            close(hw2)
        
            %removing big variables from memory
            tofsignal=nan;
            tofsignalPart_UMR=nan;
            
            tNormalTRes=h5read(dataSource,'/tofData/time');

            dataMatrix_HR=[];

        else % HR
            hrdata=readHrSpectra(dataSource);
            tNormalTRes=datenum(hrdata.time)';
        
            compoundNames = hrdata.peaks.compound;
        
            if useGPU
                dataMatrix_HR=gpuArray(hrdata.peaks.HR);
                mzs_HR=gpuArray(single(hrdata.peaks.mass));
            else
                dataMatrix_HR=hrdata.peaks.HR;
                mzs_HR=hrdata.peaks.mass;
            end
                
            IhrMassesToSkip = ismember(compoundNames,hrMassesToSkip);
            dataMatrix_HR=dataMatrix_HR(:,~IhrMassesToSkip);
            mzs_HR=mzs_HR(~IhrMassesToSkip);
            compoundNames=compoundNames(~IhrMassesToSkip);
        
            hrdata.peaks.compound = hrdata.peaks.compound(~IhrMassesToSkip);
            hrdata.peaks.HR = hrdata.peaks.HR(:,~IhrMassesToSkip);
            hrdata.peaks.mass = hrdata.peaks.mass(~IhrMassesToSkip);
            hrdata.peaks.UMR = hrdata.peaks.UMR(~IhrMassesToSkip);
            hrdata.peaks.elements = hrdata.peaks.elements(~IhrMassesToSkip,:);
        
        end

    else % for direct data import
        if useGPU
            mzs_bin = gpuArray(dataSource.mzs_bin);
            mzs_UMR = gpuArray(dataSource.mzs_UMR);
            mzs_HR = gpuArray(dataSource.mzs_HR);
            dataMatrix_bin = gpuArray(single(dataSource.dataMatrix_bin));
            dataMatrix_UMR = gpuArray(single(dataSource.dataMatrix_UMR));
            dataMatrix_HR = gpuArray(single(dataSource.dataMatrix_HR));
        else
            mzs_bin = dataSource.mzs_bin;
            mzs_UMR = dataSource.mzs_UMR;
            mzs_HR = dataSource.mzs_HR;
            dataMatrix_bin = dataSource.dataMatrix_bin;
            dataMatrix_UMR = dataSource.dataMatrix_UMR;
            dataMatrix_HR = dataSource.dataMatrix_HR;
        end

        compoundNames=dataSource.compoundNames;

        IhrMassesToSkip = ismember(compoundNames,hrMassesToSkip);
        dataMatrix_HR=dataMatrix_HR(:,~IhrMassesToSkip);
        mzs_HR=mzs_HR(~IhrMassesToSkip);
        compoundNames=compoundNames(~IhrMassesToSkip);

        tNormalTRes = dataSource.t;
    end

    tLowerTRes=tNormalTRes(1:timeResolutionLoweringFactor:end);
    
    timerange=1:size(tNormalTRes,2);
    timerangeLowerTRes=1:size(tLowerTRes,2);
    

    %for normal time res
        timeIsOk=false(size(timerange));
        for ii=1:size(goodTimes,1)
            I=tNormalTRes>=datenum(goodTimes(ii,1))-aveTimeMins/60/24&tNormalTRes<datenum(goodTimes(ii,2)); 
            timeIsOk(I)=true;
        end

        % calculating the error matrices using the defined errorParameters
        % with the sqrt(errorParameter1^2*dataMatrix+errorParameter2^2) formula
        errorMatrix_bin=single(sqrt(errorParameter1_bin^2*max(dataMatrix_bin,0) + errorParameter2_bin^2)); 
        errorMatrix_UMR=single(sqrt(errorParameter1_UMR^2*max(dataMatrix_UMR,0) + errorParameter2_UMR^2));
        errorMatrix_HR=single(sqrt(errorParameter1_HR^2*max(dataMatrix_HR,0) + errorParameter2_HR^2));


        timerange=timerange(timeIsOk); 

        if ~isempty(dataMatrix_bin)
            dataMatrix_bin=dataMatrix_bin(timerange,:);
            errorMatrix_bin=errorMatrix_bin(timerange,:);
        end

        if ~isempty(dataMatrix_UMR)
            dataMatrix_UMR=dataMatrix_UMR(timerange,:);
            errorMatrix_UMR=errorMatrix_UMR(timerange,:);
        end

        tNormalTRes=tNormalTRes(timerange);

        if strcmp(pmfType,'HR')
            dataMatrix_HR=dataMatrix_HR(timerange,:);
            errorMatrix_HR=errorMatrix_HR(timerange,:);
        end
        
        
    %for lower time res
        timeIsOk=false(size(timerangeLowerTRes));
        for ii=1:size(goodTimes,1)
            I=tLowerTRes>=datenum(goodTimes(ii,1))-(aveTimeMins*timeResolutionLoweringFactor)/60/24&tLowerTRes<datenum(goodTimes(ii,2)); 
            timeIsOk(I)=true;
        end

       
        timerangeLowerTRes=timerangeLowerTRes(timeIsOk);
        if ~isempty(dataMatrix_bin)
            dataMatrixLowerTRes=smoothi(dataMatrix_bin,timeResolutionLoweringFactor);
           
            % calculating the error matrix using the defined errorParameters
            % with the sqrt(errorParameter1^2*dataMatrix+errorParameter2^2) formula
            errorMatrixLowerTRes=single(sqrt(errorParameter1_bin^2*max(dataMatrixLowerTRes,0)/sqrt(timeResolutionLoweringFactor) + errorParameter2_bin^2)); 
    
            dataMatrixLowerTRes=dataMatrixLowerTRes(timerangeLowerTRes,:);
            tLowerTRes=tLowerTRes(timerangeLowerTRes);
            errorMatrixLowerTRes=errorMatrixLowerTRes(timerangeLowerTRes,:);
        end

    Smin=1e-6; % the minimum error possible, used to prevent zeros

    clc

end % initialization of the variables ends

% downweighting for binPMF

% for normal time resolution
    Iok=true(size(dataMatrix_bin,1),1);
    X_meas_bin=dataMatrix_bin;

    S_bin=errorMatrix_bin;
    S_bin=S_bin(Iok,:);
    S_bin(S_bin<Smin)=Smin;

    X_meas_bin=X_meas_bin(Iok,:);
    X_meas_bin(X_meas_bin<0|isnan(X_meas_bin))=0;
    SNR_bin=single(X_meas_bin./S_bin);

    if downWeightingBinOn
        S_bin(SNR_bin<2)=S_bin(SNR_bin<2)*2; %downweighting with the factor of 2
        S_bin(SNR_bin<0.2)=S_bin(SNR_bin<0.2)*5; %downweighting with the factor of 10 (because the previous line already downweights with the factor of 2)
        nDownWeigted_bin=single(zeros(size(X_meas_bin)));
        nDownWeigted_bin(SNR_bin<2)=0.5;
        nDownWeigted_bin(SNR_bin<0.2)=0.9;

        nDownWeigted_bin=sum(sum(nDownWeigted_bin));
    else
        nDownWeigted_bin=0;
    end


    % for lower time resolution
        if ~isempty(dataMatrix_bin)
            Iok=true(size(dataMatrixLowerTRes,1),1);
        
            S_binLowerTRes=errorMatrixLowerTRes;
            S_binLowerTRes=S_binLowerTRes(Iok,:);
            S_binLowerTRes(S_binLowerTRes<Smin)=Smin;
        
            X_meas_binLowerTRes=dataMatrixLowerTRes;
            X_meas_binLowerTRes=X_meas_binLowerTRes(Iok,:);
            X_meas_binLowerTRes(X_meas_binLowerTRes<0|isnan(X_meas_binLowerTRes))=0;
        
            SNR_binLowerTRes=single(X_meas_binLowerTRes./S_binLowerTRes);
            S_binLowerTRes(SNR_binLowerTRes<2)=S_binLowerTRes(SNR_binLowerTRes<2)*2; %downweighting with the factor of 2
            S_binLowerTRes(SNR_binLowerTRes<0.2)=S_binLowerTRes(SNR_binLowerTRes<0.2)*5; %downweighting with the factor of 10 (because the previous line already downweights with the factor of 2)
            nDownWeigted_binLowerTRes=single(zeros(size(X_meas_binLowerTRes)));
            nDownWeigted_binLowerTRes(SNR_binLowerTRes<2)=0.5;
            nDownWeigted_binLowerTRes(SNR_binLowerTRes<0.2)=0.9;
            nDownWeigted_binLowerTRes=sum(sum(nDownWeigted_binLowerTRes));
        end
    
% downweighting for UMR PMF

% for normal time resolution
    Iok=true(size(dataMatrix_UMR,1),1);

    S_UMR=errorMatrix_UMR;
    S_UMR=S_UMR(Iok,:);
    S_UMR(S_UMR<Smin)=Smin;

    X_meas_UMR=dataMatrix_UMR;
    X_meas_UMR=X_meas_UMR(Iok,:);
    X_meas_UMR(X_meas_UMR<0|isnan(X_meas_UMR))=0;

    SNR_UMR=single(X_meas_UMR./S_UMR);
    S_UMR(SNR_UMR<2)=S_UMR(SNR_UMR<2)*2; %downweighting with the factor of 2
    S_UMR(SNR_UMR<0.2)=S_UMR(SNR_UMR<0.2)*5; %downweighting with the factor of 10 (because the previous line already downweights with the factor of 2)
    nDownWeigted_UMR=single(zeros(size(X_meas_UMR)));
    nDownWeigted_UMR(SNR_UMR<2)=0.5;
    nDownWeigted_UMR(SNR_UMR<0.2)=0.9;
    nDownWeigted_UMR=sum(sum(nDownWeigted_UMR));

% for HR PMF

%for normal time resolution
    if strcmp(pmfType,'HR')
        Iok=true(size(dataMatrix_HR,1),1);
    
        S_HR=errorMatrix_HR;
        S_HR=S_HR(Iok,:);
        S_HR(S_HR<Smin)=Smin;
    
        X_meas_HR=dataMatrix_HR;
        X_meas_HR=X_meas_HR(Iok,:);
        X_meas_HR(X_meas_HR<0|isnan(X_meas_HR))=0;
    
        SNR_HR=single(X_meas_HR./S_HR);
        S_HR(SNR_HR<2)=S_HR(SNR_HR<2)*2; %downweighting with the factor of 2
        S_HR(SNR_HR<0.2)=S_HR(SNR_HR<0.2)*5; %downweighting with the factor of 10 (because the previous line already downweights with the factor of 2)
        nDownWeigted_HR=single(zeros(size(X_meas_HR)));
        nDownWeigted_HR(SNR_HR<2)=0.5;
        nDownWeigted_HR(SNR_HR<0.2)=0.9;
        nDownWeigted_HR=sum(sum(nDownWeigted_HR));
    end

close all

if strcmp(pmfType,'UMR')
    forceUMR = 1; %if 1, code calculates in UMR
else
    forceUMR = 0; %if 1, code calculates in UMR
end

% initializing vectors to track Q and unexplained variation
QQexp_min_vec=[];
QQexp_robust_min_vec=[];
QQexp_true_min_vec=[];
unExplainedVariation_vec=[];
p_vec=[];

plotpmfstats

for p=pRange
    
    QQexp_min=inf;
    for iRepetition=1:repetitions
            
            waitbartext={['Running PMF with p = ' int2str(p)],'',['Selected p range: ' int2str(pRange)],['Repetitions per p: ' int2str(repetitions)]};
            if ishandle(hw)
                waitbar((iRepetition-1)/repetitions,hw,waitbartext);
            else
                hw=waitbar((iRepetition-1)/repetitions,waitbartext);
                hw.Position=[1100 700 270 100];
                hw.Name='PMF running';
                drawnow
            end
        
        if forceUMR % UMR
        
            X_meas=X_meas_UMR; 
            X_meas_nonConstrained=X_meas_UMR; 
            S=S_UMR;
            nDownWeigted=nDownWeigted_UMR;
            t=tNormalTRes;
        
            if nConstrainedFactors>0
                X_meas=max(0,X_meas_UMR-constrainedTS*constrainedMS);
            end
        
            mzs=mzs_UMR;
            binWidth=binWidth_UMR;
            if isnan(binWidth)
                binWidth=1;
            end
        
            umrQQexpTrackedWithBinPMF=0;
        else
        
            if ~strcmp(pmfType,'HR') % bin
                X_meas=X_meas_bin;
                X_meas_nonConstrained=X_meas_bin;
                S=S_bin;
                nDownWeigted=nDownWeigted_bin;
                mzs=mzs_bin;
                binWidth=binWidth_bin;

                if nConstrainedFactors>0
                    X_meas=max(0,X_meas_bin-constrainedTS*constrainedMS);
                end

            else % HR
                X_meas=X_meas_HR;
                X_meas_nonConstrained=X_meas_HR;
                S=S_HR;
                nDownWeigted=nDownWeigted_HR;
                mzs=mzs_HR;
                binWidth=0;

                if nConstrainedFactors>0
                    X_meas=max(0,X_meas_HR-constrainedTS*constrainedMS);
                end

            end

            t=tNormalTRes;
        
        end
        
        TS_UMR=[];
        MS_UMR=[];
        if strcmp(pmfAlg,'PSPMF')
             [TS,MS,TS_UMR,MS_UMR]=pspmfinternal(pspmfAlgs,pspmfRepeats,...
                 {X_meas_UMR,X_meas_binLowerTRes,X_meas_bin},p,...
                 {S_UMR,S_binLowerTRes,S_bin},...
                 {nDownWeigted_UMR,nDownWeigted_binLowerTRes,nDownWeigted_bin});
        
        else
            if nConstrainedFactors==0
                [TS,MS]=pmfinternal_algs(pmfAlg,X_meas,p,S,nDownWeigted);
            else
                [TS,MS]=pmfinternal_algs_con(pmfAlg,X_meas,p,S,nDownWeigted,X_meas_nonConstrained,nConstrainedFactors,involveConstrainedFactorsInQCalc,constrainedTS,constrainedMS);
            end
        end
        
        if ishandle(1)
            
        else
            break
        end
                
        E=X_meas-TS*MS; % residual
        
        if robustModeOn
            robustnessNumber=16; % to use the robust method with abs > 4 forced to 4
        else
            robustnessNumber=inf;
        end
        
        % calculate Q/Qexp values in a robust way and a true way
        if ~isempty(TS_UMR)&&~isempty(MS_UMR) % for PSPMF operation
            E_UMR=X_meas_UMR-TS_UMR*MS_UMR;
            QQexp_UMR=sum(sum(min((E_UMR./S_UMR).^2,robustnessNumber)))/(numel(X_meas_UMR)-nDownWeigted_UMR-numel(TS_UMR)-numel(MS_UMR));
            QQexp_robust_UMR=sum(sum(min((E_UMR./S_UMR).^2,16)))/(numel(X_meas_UMR)-nDownWeigted_UMR-numel(TS_UMR)-numel(MS_UMR));
            QQexp_true_UMR=sum(sum((E_UMR./S_UMR).^2))/(numel(X_meas_UMR)-nDownWeigted_UMR-numel(TS_UMR)-numel(MS_UMR));
        
        else
            E_UMR=[];
            QQexp_UMR=[];
            QQexp_robust_UMR=[];
            QQexp_true_UMR=[];
        end
        
        QQexp=sum(sum(min((E./S).^2,robustnessNumber)))/(numel(X_meas)-nDownWeigted-numel(TS)-numel(MS));
        QQexp_robust=sum(sum(min((E./S).^2,16)))/(numel(X_meas)-nDownWeigted-numel(TS)-numel(MS));
        QQexp_true=sum(sum((E./S).^2))/(numel(X_meas)-nDownWeigted-numel(TS)-numel(MS));
        
        if nConstrainedFactors>0
            E=X_meas_nonConstrained-[TS constrainedTS] * [MS; constrainedMS];
        
            QQexp=sum(sum(min((E./S).^2,robustnessNumber)))/(numel(X_meas_nonConstrained)-nDownWeigted-(p+nConstrainedFactors)*(size(X_meas_nonConstrained,1)+size(X_meas_nonConstrained,2)));
            QQexp_robust=sum(sum(min((E./S).^2,16)))/(numel(X_meas_nonConstrained)-nDownWeigted-(p+nConstrainedFactors)*(size(X_meas_nonConstrained,1)+size(X_meas_nonConstrained,2)));
            QQexp_true=sum(sum((E./S).^2))/(numel(X_meas_nonConstrained)-nDownWeigted-(p+nConstrainedFactors)*(size(X_meas_nonConstrained,1)+size(X_meas_nonConstrained,2)));
        end
        
        if umrQQexpTrackedWithBinPMF
            QQexp=QQexp_UMR;
            QQexp_robust=QQexp_robust_UMR;
            QQexp_true=QQexp_true_UMR;
        end
           
        % when a better solution is found            
        if QQexp<QQexp_min
            QQexp_min=QQexp;
        
            QQexp_robust_min=QQexp_robust;
            QQexp_true_min=QQexp_true;
                    
            [~,I]=sort(sum(MS.*(mzs+binWidth/2),2)./sum(MS,2)); % sort factors according to their average masses
            TS=TS(:,I);
            MS=MS(I,:);
            
            % normalize factors so that their MS sums to 1
            MSsum=sum(MS,2);
            MS_final=MS./MSsum;
            TS_final=TS.*MSsum';
            
            QQexp_final=QQexp;
            E_final=E;
        
            % plot the currently best solution in real time
            plotpmfresults

            % store the currently best solution as variables
            if forceUMR % UMR

                eval(['UMRPMFRes' int2str(p) 'Factors.MS=MS_final;'])
                eval(['UMRPMFRes' int2str(p) 'Factors.TS=TS_final;'])
                eval(['UMRPMFRes' int2str(p) 'Factors.t=t;'])
                eval(['UMRPMFRes' int2str(p) 'Factors.mzs=mzs;'])
            else

                if ~strcmp(pmfType,'HR') % bin
        
                    eval(['binPMFRes' int2str(p) 'Factors.MS=MS_final;'])
                    eval(['binPMFRes' int2str(p) 'Factors.TS=TS_final;'])
                    eval(['binPMFRes' int2str(p) 'Factors.t=t;'])
                    eval(['binPMFRes' int2str(p) 'Factors.mzs=mzs;'])

                else % HR

                    eval(['HRPMFRes' int2str(p) 'Factors.MS=MS_final;'])
                    eval(['HRPMFRes' int2str(p) 'Factors.TS=TS_final;'])
                    eval(['HRPMFRes' int2str(p) 'Factors.t=t;'])
                    eval(['HRPMFRes' int2str(p) 'Factors.mzs=mzs;'])
                    eval(['HRPMFRes' int2str(p) 'Factors.compoundNames=compoundNames;'])

                end
            end
        end
    end
    
    close(hw)

    % update Q and unexplained variation data for tracking
    QQexp_min_vec=[QQexp_min_vec QQexp_min];
    QQexp_robust_min_vec=[QQexp_robust_min_vec QQexp_robust_min];
    QQexp_true_min_vec=[QQexp_true_min_vec QQexp_true_min];
    unExplainedVariation_vec =[unExplainedVariation_vec unExplainedVariation];
    p_vec=[p_vec p];
    
    % plot the Q/Qexp vs p plot in real time
    if ishandle(1)
        plotpmfstats
    else
        plotpmfstats
        break
    end

end

if ishandle(hw)
    close(hw)
end

figure(1)
subplot(211)
title('')