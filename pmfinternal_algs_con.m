function [TS,MS]=pmfinternal_algs_con(pmfalg,X,p,S,nDownWeigted,X_full,nConstrainedFactors,involveConstrainedFactorsInQCalc,constrainedTS,constrainedMS)


if strcmp(pmfalg,'ALS')||strcmp(pmfalg,'MULT')
    [TS,MS]=pmfinternal_con(X,p,S,nDownWeigted,1,pmfalg,[],[],X_full,nConstrainedFactors,involveConstrainedFactorsInQCalc,constrainedTS,constrainedMS);

else
    [TS,MS]=pmfinternal_con(X,p,S,nDownWeigted,str2double(pmfalg(4)),'ALS',[],[],X_full,nConstrainedFactors,involveConstrainedFactorsInQCalc,constrainedTS,constrainedMS);
    [TS,MS]=pmfinternal_con(X,p,S,nDownWeigted,1,'MULT',TS,MS,X_full,nConstrainedFactors,involveConstrainedFactorsInQCalc,constrainedTS,constrainedMS);

end
