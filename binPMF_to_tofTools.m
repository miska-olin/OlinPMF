%% Define the new data file name for tofTools

tofToolsDataFileNameFromBinPMF='.\data\paramEvalFromBinPMF.h5';

%% processing the already loaded data
parsnew=[mean(pars(1,pars(1,:)>0));mean(pars(2,pars(2,:)<-1e4))];
tofs=funcArgs.inverseEQ(mzs_bin+binWidth_bin/2,parsnew);
signalLength=h5info(dataSource,'/tofData/signal').Dataspace.Size(1);

clear tofsignalmod
for i=1:p
    tofsignalmod(:,i)=interp1(tofs,MS_final(i,:),1:signalLength);
end

%% copying the raw tofTools file as the new file
system(['copy ' dataSource ' ' tofToolsDataFileNameFromBinPMF]);

fid = H5F.open(tofToolsDataFileNameFromBinPMF,'H5F_ACC_RDWR','H5P_DEFAULT');
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

tof_write_HDF(tofToolsDataFileNameFromBinPMF,'/tofData/SI',repmat(mean(SIorig,'omitnan'),p,1),'memType','H5T_NATIVE_FLOAT','chunk',[p 1]);
tof_write_HDF(tofToolsDataFileNameFromBinPMF,'/tofData/baseline',zeros(p,signalLength),'memType','H5T_NATIVE_FLOAT','chunk',[1 signalLength]);
tof_write_HDF(tofToolsDataFileNameFromBinPMF,'/tofData/coadds', median(coadds(coadds>0))*ones(p,1),'memType','H5T_NATIVE_FLOAT','chunk',[p 1]);
tof_write_HDF(tofToolsDataFileNameFromBinPMF,'/tofData/signal',tofsignalmod','memType','H5T_NATIVE_FLOAT');
tof_write_HDF(tofToolsDataFileNameFromBinPMF,'/tofData/time',(timeorig(1)+(0:(p-1))*aveTimeMins/60/24)','memType','H5T_NATIVE_DOUBLE','chunk',[p 1]);

tof_write_HDF(tofToolsDataFileNameFromBinPMF,'/MassCalib/a_log',(repmat(mean(alogorig,'omitnan'),1,p))','memType','H5T_NATIVE_FLOAT','chunk',[p 1]);
tof_write_HDF(tofToolsDataFileNameFromBinPMF,'/MassCalib/b_log',(repmat(mean(blogorig,'omitnan'),1,p))','memType','H5T_NATIVE_FLOAT','chunk',[p 1]);
tof_write_HDF(tofToolsDataFileNameFromBinPMF,'/MassCalib/masses',(repmat(massesorig(:,1),1,p))','memType','H5T_NATIVE_FLOAT','chunk',[p nMassCalibrants]);
tof_write_HDF(tofToolsDataFileNameFromBinPMF,'/MassCalib/pars',(repmat(parsnew,1,p))','memType','H5T_NATIVE_FLOAT','chunk',[p 2]);
tof_write_HDF(tofToolsDataFileNameFromBinPMF,'/MassCalib/ppm',(repmat(1,1,p))','memType','H5T_NATIVE_FLOAT','chunk',[p 1]);
tof_write_HDF(tofToolsDataFileNameFromBinPMF,'/MassCalib/ps_log',(repmat(2,1,p))','memType','H5T_NATIVE_FLOAT','chunk',[p 1]);
tof_write_HDF(tofToolsDataFileNameFromBinPMF,'/MassCalib/sn',(repmat([mean(snorig(1,snorig(1,:)>0));mean(snorig(2,snorig(2,:)>0));mean(snorig(3,snorig(3,:)>0))],1,p))','memType','H5T_NATIVE_FLOAT','chunk',[p nMassCalibrants]);



tof_write_attributes(tofToolsDataFileNameFromBinPMF,{'Noise reduction time step [min]',0},'/tofData');
tof_write_attributes(tofToolsDataFileNameFromBinPMF,{'Averaging time step [min]',aveTimeMins},'/tofData');
tof_write_attributes(tofToolsDataFileNameFromBinPMF,{'Signal Unit','mV/ext'},'/tofData');
tof_write_attributes(tofToolsDataFileNameFromBinPMF,{'Start Time',h5readatt(dataSource,'/tofData','Start Time')},'/tofData');
tof_write_attributes(tofToolsDataFileNameFromBinPMF,{'Stop Time',h5readatt(dataSource,'/tofData','Start Time')+p*aveTimeMins/60/24},'/tofData');
tof_write_attributes(tofToolsDataFileNameFromBinPMF,{'Tof Frequency',h5readatt(dataSource,'/tofData','Tof Frequency')},'/tofData');
tof_write_attributes(tofToolsDataFileNameFromBinPMF,{'Sampling Interval (s)',h5readatt(dataSource,'/tofData','Sampling Interval (s)')},'/tofData');
tof_write_attributes(tofToolsDataFileNameFromBinPMF,{'Polarity','positive'},'/tofData');
tof_write_attributes(tofToolsDataFileNameFromBinPMF,{'baselineCorrection',1},'/tofData');
tof_write_attributes(tofToolsDataFileNameFromBinPMF,{'hasIMS',0},'/tofData');
tof_write_attributes(tofToolsDataFileNameFromBinPMF,{'diffSpec',0},'/tofData');
tof_write_attributes(tofToolsDataFileNameFromBinPMF,{'DC_card','Tuatara'},'/tofData');
tof_write_attributes(tofToolsDataFileNameFromBinPMF,{'NdigoProcessingType',-1},'/tofData');
tof_write_attributes(tofToolsDataFileNameFromBinPMF,{'UsedProcessingType','ADC'},'/tofData');

tof_write_attributes(tofToolsDataFileNameFromBinPMF,{'names',cellfun(@(x)x(1:end-1),cellstr(h5readatt(dataSource,'/MassCalib','names')),'UniformOutput',false)},'/MassCalib');
tof_write_attributes(tofToolsDataFileNameFromBinPMF,{'masses',{h5readatt(dataSource,'/MassCalib','masses')}},'/MassCalib');
tof_write_attributes(tofToolsDataFileNameFromBinPMF,{'equation','tof_fit_tof_mz'},'/MassCalib');
tof_write_attributes(tofToolsDataFileNameFromBinPMF,{'mxNrMasses',nMassCalibrants},'/MassCalib');
tof_write_attributes(tofToolsDataFileNameFromBinPMF,{'timeRes',h5readatt(dataSource,'/MassCalib','timeRes')},'/MassCalib');
tof_write_attributes(tofToolsDataFileNameFromBinPMF,{'polarity','positive'},'/MassCalib');
tof_write_attributes(tofToolsDataFileNameFromBinPMF,{'doPlot',0},'/MassCalib');
tof_write_attributes(tofToolsDataFileNameFromBinPMF,{'mode','moving'},'/MassCalib');
tof_write_attributes(tofToolsDataFileNameFromBinPMF,{'a',parsnew(1)},'/MassCalib');
tof_write_attributes(tofToolsDataFileNameFromBinPMF,{'b',parsnew(2)},'/MassCalib');
tof_write_attributes(tofToolsDataFileNameFromBinPMF,{'peakShape',1},'/MassCalib');
tof_write_attributes(tofToolsDataFileNameFromBinPMF,{'useLoggedParameters',0},'/MassCalib');
tof_write_attributes(tofToolsDataFileNameFromBinPMF,{'sparse',0},'/MassCalib');
tof_write_attributes(tofToolsDataFileNameFromBinPMF,{'peakShapeId',char(h5readatt(dataSource,'/MassCalib','peakShapeId'))},'/MassCalib');

disp("binPMF solution has been succesfully written into " + tofToolsDataFileNameFromBinPMF + " file.")
%% make the file smaller

if isunix
    exitFlag=system("h5repack " + strrep(tofToolsDataFileNameFromBinPMF,filesep,'/') + " " + strrep(tofToolsDataFileNameFromBinPMF(1:end-3),filesep,'/') + "Small.h5");
    if exitFlag==0 % success
        system("rm " + tofToolsDataFileNameFromBinPMF);
        disp("The file has been made smaller succesfully by removing unnecessary data fields included in the original datafile.")
        disp("The file name has been changed to end with Small.h5")
    else
        disp("Failed to make the file smaller due to missing h5repack command.")
        disp("The file is unnecessary large because it includes also data fields included in the original datafile.")
    end
else % windows
    exitFlag=system("wsl h5repack " + strrep(tofToolsDataFileNameFromBinPMF,filesep,'/') + " " + strrep(tofToolsDataFileNameFromBinPMF(1:end-3),filesep,'/') + "Small.h5");
    if exitFlag==0 % success
        system("del " + tofToolsDataFileNameFromBinPMF);
        disp("The file has been made smaller succesfully by removing unnecessary data fields included in the original datafile.")
        disp("The file name has been changed to end with Small.h5")
    else
        disp("Failed to make the file smaller due to missing wsl and/or h5repack commands.")
        disp("The file is unnecessary large because it includes also data fields included in the original datafile.")
    end
end
