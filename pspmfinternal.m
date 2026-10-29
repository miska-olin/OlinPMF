function [TS,MS,TS_UMR,MS_UMR]=pspmfinternal(pspmfAlgs,pspmfRepeats,X_meas_s,p,Ss,nDownWeigteds)


% UMR normal time resolution run
    i=1;
    disp(['PSPMF: UMR normal time resolution run with ' pspmfAlgs{i} ' algorithm(s) having ' int2str(pspmfRepeats{i}) ' repeats'])
    
    nAlgs=length(pspmfAlgs{i});
    
    iAlg=1;
    if strcmp(pspmfAlgs{i}(iAlg),'A')
        [TS,MS]=pmfinternal(X_meas_s{i},p,Ss{i},nDownWeigteds{i},pspmfRepeats{i}(iAlg),'ALS',[],[]);

    elseif strcmp(pspmfAlgs{i}(iAlg),'M')
        [TS,MS]=pmfinternal(X_meas_s{i},p,Ss{i},nDownWeigteds{i},pspmfRepeats{i}(iAlg),'MULT',[],[]);
    end
    
    if nAlgs==2
        iAlg=2;
        if strcmp(pspmfAlgs{i}(iAlg),'A')
            [TS,MS]=pmfinternal(X_meas_s{i},p,Ss{i},nDownWeigteds{i},pspmfRepeats{i}(iAlg),'ALS',TS,MS);
        elseif strcmp(pspmfAlgs{i}(iAlg),'M')
            [TS,MS]=pmfinternal(X_meas_s{i},p,Ss{i},nDownWeigteds{i},pspmfRepeats{i}(iAlg),'MULT',TS,MS);
        end
    end

    TSi=TS;
    TS_UMR=TS;
    MS_UMR=MS;

% bin lower time resolution run
    i=2;
    disp(['PSPMF: bin lower time resolution run with ' pspmfAlgs{i} ' algorithm(s) having ' int2str(pspmfRepeats{i}) ' repeats'])
    
    nAlgs=length(pspmfAlgs{i});
    
    
    iAlg=1;
    if pspmfRepeats{i}(iAlg)>0
        if strcmp(pspmfAlgs{i}(iAlg),'A')
            [TS,MS]=pmfinternal(X_meas_s{i},p,Ss{i},nDownWeigteds{i},pspmfRepeats{i}(iAlg),'ALS',[],[]);
        elseif strcmp(pspmfAlgs{i}(iAlg),'M')
            [TS,MS]=pmfinternal(X_meas_s{i},p,Ss{i},nDownWeigteds{i},pspmfRepeats{i}(iAlg),'MULT',[],[]);
        end
    else
        TS=[];
        MS=[];
    end
    
    if nAlgs==2
        iAlg=2;
        if strcmp(pspmfAlgs{i}(iAlg),'A')
            [TS,MS]=pmfinternal(X_meas_s{i},p,Ss{i},nDownWeigteds{i},pspmfRepeats{i}(iAlg),'ALS',TS,MS);
        elseif strcmp(pspmfAlgs{i}(iAlg),'M')
            [TS,MS]=pmfinternal(X_meas_s{i},p,Ss{i},nDownWeigteds{i},pspmfRepeats{i}(iAlg),'MULT',TS,MS);
        end
    end

    MSi=MS;

% bin normal time resolution run
    i=3;
    disp(['PSPMF: bin normal time resolution run (final) with ' pspmfAlgs{i} ' algorithm(s) having ' int2str(pspmfRepeats{i}) ' repeats'])
    
    nAlgs=length(pspmfAlgs{i});
    
    iAlg=1;
    if strcmp(pspmfAlgs{i}(iAlg),'A')
        [TS,MS]=pmfinternal(X_meas_s{i},p,Ss{i},nDownWeigteds{i},pspmfRepeats{i}(iAlg),'ALS',TSi,MSi);
    elseif strcmp(pspmfAlgs{i}(iAlg),'M')
        [TS,MS]=pmfinternal(X_meas_s{i},p,Ss{i},nDownWeigteds{i},pspmfRepeats{i}(iAlg),'MULT',TSi,MSi);
    end
    
    if nAlgs==2
        iAlg=2;
        if strcmp(pspmfAlgs{i}(iAlg),'A')
            [TS,MS]=pmfinternal(X_meas_s{i},p,Ss{i},nDownWeigteds{i},pspmfRepeats{i}(iAlg),'ALS',TS,MS);
        elseif strcmp(pspmfAlgs{i}(iAlg),'M')
            [TS,MS]=pmfinternal(X_meas_s{i},p,Ss{i},nDownWeigteds{i},pspmfRepeats{i}(iAlg),'MULT',TS,MS);
        end
    end

