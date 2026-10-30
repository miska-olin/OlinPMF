%% Define the new RR-data file name for tofTools and a temporary MAT-file

tofToolsDataFileNameFromRRPMF='.\data\paramEvalFromRRPMF.h5';
RRDataMatFileName = '.\data\RRdata.mat';

%% Define parameters for RR-PMF 

initializeRRVariables = 0; % Set to 1 to make the temporary MAT-file. 
% Note that this step may take a long time, depending on your data size.
% If the parameters below will change, the initialization has to be re-run.

mzRange = [40:45 47:54 58:72 74:400]; % what masses are included (always listed as UMR-type masses)
binWidth_RR = 0.001; % width of a bin in RR-PMF (in Th)
lowerMassDefect = 0.25;  % UMR minus this value (in Th) is the lower limit of a mass
upperMassDefect = 0.35; % UMR plus this value (in Th) is the upper limit of a mass

% Note that TS_final matrix must contain the UMR-PMF solution that you want
% RR-PMF to apply the Half-ALS algorithm. If TS_final matrix is not the
% intended UMR-PMF solution, you can load the solution there by command,
% e.g., loadpmfresults(5,'UMR')

%% making the MAT-file if requested

mzs_RR=[];
lowerBins=floor(lowerMassDefect/binWidth_RR);
upperBins=floor(upperMassDefect/binWidth_RR);

for i=1:length(mzRange)
    mzs_RR=[mzs_RR mzRange(i)];

    for j=1:lowerBins
        mzs_RR=[mzs_RR mzRange(i)-j*binWidth_RR];
    end

    for j=1:upperBins-1
        mzs_RR=[mzs_RR mzRange(i)+j*binWidth_RR];
    end
end
mzs_RR=sort(mzs_RR);


convcoef=tof_get_convcoef(dataSource);
[pars,eq] = tof_get_mz_calib_pars(dataSource);
funcArgs = eq();

if initializeRRVariables
    disp('Writing tofTools RR-data to a MAT-file')

    RRData=matfile(RRDataMatFileName,'Writable',true);
    RRData.X_meas_RR=single([]);
    
    linestoread=1000;
    
    timevalidities=false(h5info(dataSource,'/tofData/signal').Dataspace.Size(2),1);
    timevalidities(timerange)=true;
    
    datasplitting1=1:linestoread:h5info(dataSource,'/tofData/signal').Dataspace.Size(2);
    datasplitting2=datasplitting1+linestoread-1;
    datasplitting2(end)=h5info(dataSource,'/tofData/signal').Dataspace.Size(2);
    toflength=h5info(dataSource,'/tofData/signal').Dataspace.Size(1);
    
    for iPart=1:length(datasplitting1)
        progressbar((iPart-1)/length(datasplitting1));
    
        if iPart==length(datasplitting1)
            linestoread=datasplitting2(end)-datasplitting1(end)+1;
        end
        
        tofsignal=h5read(dataSource,'/tofData/signal',[1 datasplitting1(iPart)],[toflength linestoread]); 
        tofsignalPart=zeros(size(tofsignal,2),length(mzs_RR));
    
        convcoefPart=convcoef(datasplitting1(iPart):datasplitting2(iPart));
        parsPart=pars(:,datasplitting1(iPart):datasplitting2(iPart));
        timevaliditiesPart=timevalidities(datasplitting1(iPart):datasplitting2(iPart));
    
        for i=1:size(tofsignalPart,2) %masses
            for j=1:size(tofsignalPart,1) %times
                startI=funcArgs.inverseEQ(mzs_RR(i),parsPart(:,j));
                endI=funcArgs.inverseEQ(mzs_RR(i)+binWidth_RR,parsPart(:,j));
                
                startI=startI+0.5;
                endI=endI+0.5;
    
                if startI>0&&endI>0
                    tofsignalPart(j,i)=(sum(tofsignal(floor(startI):(ceil(endI)-1),j))+...
                        (floor(startI)-startI)*tofsignal(floor(startI),j)+...
                        (endI-ceil(endI))*tofsignal(ceil(endI)-1,j))*convcoefPart(j);
                end
            end
        end
    
        RRData.X_meas_RR=[RRData.X_meas_RR; tofsignalPart(timevaliditiesPart,:)];
    end
    progressbar(1);
                                                                                                
else
    RRData=matfile(RRDataMatFileName,'Writable',false);
end

%% RR-PMF step
disp('Applying Half-ALS algorithm on the data to get RR-PMF factors')
MS_final_RR=max(0,TS_final\RRData.X_meas_RR); 

%% processing the already loaded data
parsnew=[mean(pars(1,pars(1,:)>0));mean(pars(2,pars(2,:)<-1e4))];
tofs=funcArgs.inverseEQ(mzs_RR+binWidth_RR/2,parsnew);
signalLength=h5info(dataSource,'/tofData/signal').Dataspace.Size(1);

clear tofsignalmod
for i=1:p
    tofsignalmod(:,i)=interp1(tofs,MS_final_RR(i,:),1:signalLength);
end

%% copying the raw tofTools file as the new file
system(['copy ' dataSource ' ' tofToolsDataFileNameFromRRPMF]);

fid = H5F.open(tofToolsDataFileNameFromRRPMF,'H5F_ACC_RDWR','H5P_DEFAULT');
plist = 'H5P_DEFAULT';
H5L.delete(fid,'tofData',plist);
H5L.delete(fid,'MassCalib',plist);

H5F.close(fid);

%% writing data to the new file
SIorig=h5read(dataSource,'/tofData/SI');
timeorig=h5read(dataSource,'/tofData/time');
coadds=h5read(dataSource,'/tofData/coadds');

alogorig=h5read(dataSource,'/MassCalib/a_log');
blogorig=h5read(dataSource,'/MassCalib/b_log');
massesorig=h5read(dataSource,'/MassCalib/masses');
snorig=h5read(dataSource,'/MassCalib/sn');
nMassCalibrants=h5readatt(dataSource,'/MassCalib','mxNrMasses');

tof_write_HDF(tofToolsDataFileNameFromRRPMF,'/tofData/SI',repmat(mean(SIorig,'omitnan'),p,1),'memType','H5T_NATIVE_FLOAT','chunk',[p 1]);
tof_write_HDF(tofToolsDataFileNameFromRRPMF,'/tofData/baseline',zeros(p,signalLength),'memType','H5T_NATIVE_FLOAT','chunk',[1 signalLength]);
tof_write_HDF(tofToolsDataFileNameFromRRPMF,'/tofData/coadds', median(coadds(coadds>0))*ones(p,1),'memType','H5T_NATIVE_FLOAT','chunk',[p 1]);
tof_write_HDF(tofToolsDataFileNameFromRRPMF,'/tofData/signal',tofsignalmod','memType','H5T_NATIVE_FLOAT');
tof_write_HDF(tofToolsDataFileNameFromRRPMF,'/tofData/time',(timeorig(1)+(0:(p-1))*aveTimeMins/60/24)','memType','H5T_NATIVE_DOUBLE','chunk',[p 1]);

tof_write_HDF(tofToolsDataFileNameFromRRPMF,'/MassCalib/a_log',(repmat(mean(alogorig,'omitnan'),1,p))','memType','H5T_NATIVE_FLOAT','chunk',[p 1]);
tof_write_HDF(tofToolsDataFileNameFromRRPMF,'/MassCalib/b_log',(repmat(mean(blogorig,'omitnan'),1,p))','memType','H5T_NATIVE_FLOAT','chunk',[p 1]);
tof_write_HDF(tofToolsDataFileNameFromRRPMF,'/MassCalib/masses',(repmat(massesorig(:,1),1,p))','memType','H5T_NATIVE_FLOAT','chunk',[p nMassCalibrants]);
tof_write_HDF(tofToolsDataFileNameFromRRPMF,'/MassCalib/pars',(repmat(parsnew,1,p))','memType','H5T_NATIVE_FLOAT','chunk',[p 2]);
tof_write_HDF(tofToolsDataFileNameFromRRPMF,'/MassCalib/ppm',(repmat(1,1,p))','memType','H5T_NATIVE_FLOAT','chunk',[p 1]);
tof_write_HDF(tofToolsDataFileNameFromRRPMF,'/MassCalib/ps_log',(repmat(2,1,p))','memType','H5T_NATIVE_FLOAT','chunk',[p 1]);
tof_write_HDF(tofToolsDataFileNameFromRRPMF,'/MassCalib/sn',(repmat([mean(snorig(1,snorig(1,:)>0));mean(snorig(2,snorig(2,:)>0));mean(snorig(3,snorig(3,:)>0))],1,p))','memType','H5T_NATIVE_FLOAT','chunk',[p nMassCalibrants]);



tof_write_attributes(tofToolsDataFileNameFromRRPMF,{'Noise reduction time step [min]',0},'/tofData');
tof_write_attributes(tofToolsDataFileNameFromRRPMF,{'Averaging time step [min]',aveTimeMins},'/tofData');
tof_write_attributes(tofToolsDataFileNameFromRRPMF,{'Signal Unit','mV/ext'},'/tofData');
tof_write_attributes(tofToolsDataFileNameFromRRPMF,{'Start Time',h5readatt(dataSource,'/tofData','Start Time')},'/tofData');
tof_write_attributes(tofToolsDataFileNameFromRRPMF,{'Stop Time',h5readatt(dataSource,'/tofData','Start Time')+p*aveTimeMins/60/24},'/tofData');
tof_write_attributes(tofToolsDataFileNameFromRRPMF,{'Tof Frequency',h5readatt(dataSource,'/tofData','Tof Frequency')},'/tofData');
tof_write_attributes(tofToolsDataFileNameFromRRPMF,{'Sampling Interval (s)',h5readatt(dataSource,'/tofData','Sampling Interval (s)')},'/tofData');
tof_write_attributes(tofToolsDataFileNameFromRRPMF,{'Polarity','positive'},'/tofData');
tof_write_attributes(tofToolsDataFileNameFromRRPMF,{'baselineCorrection',1},'/tofData');
tof_write_attributes(tofToolsDataFileNameFromRRPMF,{'hasIMS',0},'/tofData');
tof_write_attributes(tofToolsDataFileNameFromRRPMF,{'diffSpec',0},'/tofData');
tof_write_attributes(tofToolsDataFileNameFromRRPMF,{'DC_card','Tuatara'},'/tofData');
tof_write_attributes(tofToolsDataFileNameFromRRPMF,{'NdigoProcessingType',-1},'/tofData');
tof_write_attributes(tofToolsDataFileNameFromRRPMF,{'UsedProcessingType','ADC'},'/tofData');

tof_write_attributes(tofToolsDataFileNameFromRRPMF,{'names',cellfun(@(x)x(1:end-1),cellstr(h5readatt(dataSource,'/MassCalib','names')),'UniformOutput',false)},'/MassCalib');
tof_write_attributes(tofToolsDataFileNameFromRRPMF,{'masses',{h5readatt(dataSource,'/MassCalib','masses')}},'/MassCalib');
tof_write_attributes(tofToolsDataFileNameFromRRPMF,{'equation','tof_fit_tof_mz'},'/MassCalib');
tof_write_attributes(tofToolsDataFileNameFromRRPMF,{'mxNrMasses',nMassCalibrants},'/MassCalib');
tof_write_attributes(tofToolsDataFileNameFromRRPMF,{'timeRes',h5readatt(dataSource,'/MassCalib','timeRes')},'/MassCalib');
tof_write_attributes(tofToolsDataFileNameFromRRPMF,{'polarity','positive'},'/MassCalib');
tof_write_attributes(tofToolsDataFileNameFromRRPMF,{'doPlot',0},'/MassCalib');
tof_write_attributes(tofToolsDataFileNameFromRRPMF,{'mode','moving'},'/MassCalib');
tof_write_attributes(tofToolsDataFileNameFromRRPMF,{'a',parsnew(1)},'/MassCalib');
tof_write_attributes(tofToolsDataFileNameFromRRPMF,{'b',parsnew(2)},'/MassCalib');
tof_write_attributes(tofToolsDataFileNameFromRRPMF,{'peakShape',1},'/MassCalib');
tof_write_attributes(tofToolsDataFileNameFromRRPMF,{'useLoggedParameters',0},'/MassCalib');
tof_write_attributes(tofToolsDataFileNameFromRRPMF,{'sparse',0},'/MassCalib');
tof_write_attributes(tofToolsDataFileNameFromRRPMF,{'peakShapeId',char(h5readatt(dataSource,'/MassCalib','peakShapeId'))},'/MassCalib');

disp("RR-PMF solution has been succesfully written into " + tofToolsDataFileNameFromRRPMF + " file.")
%% make the file smaller

if isunix
    exitFlag=system("h5repack " + strrep(tofToolsDataFileNameFromRRPMF,filesep,'/') + " " + strrep(tofToolsDataFileNameFromRRPMF(1:end-3),filesep,'/') + "Small.h5");
    if exitFlag==0 % success
        system("rm " + tofToolsDataFileNameFromRRPMF);
        disp("The file has been made smaller succesfully by removing unnecessary data fields included in the original datafile.")
        disp("The file name has been changed to end with Small.h5")
    else
        disp("Failed to make the file smaller due to missing h5repack command.")
        disp("The file is unnecessary large because it includes also data fields included in the original datafile.")
    end
else % windows
    exitFlag=system("wsl h5repack " + strrep(tofToolsDataFileNameFromRRPMF,filesep,'/') + " " + strrep(tofToolsDataFileNameFromRRPMF(1:end-3),filesep,'/') + "Small.h5");
    if exitFlag==0 % success
        system("del " + tofToolsDataFileNameFromRRPMF);
        disp("The file has been made smaller succesfully by removing unnecessary data fields included in the original datafile.")
        disp("The file name has been changed to end with Small.h5")
    else
        disp("Failed to make the file smaller due to missing wsl and/or h5repack commands.")
        disp("The file is unnecessary large because it includes also data fields included in the original datafile.")
    end
end
