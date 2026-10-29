function [TS,MS]=pmfinternal_algs(pmfalg,X,p,S,nDownWeigted)


if strcmp(pmfalg,'ALS')||strcmp(pmfalg,'MULT')
    [TS,MS]=pmfinternal(X,p,S,nDownWeigted,1,pmfalg,[],[]);

else
    [TS,MS]=pmfinternal(X,p,S,nDownWeigted,str2double(pmfalg(4)),'ALS',[],[]);
    [TS,MS]=pmfinternal(X,p,S,nDownWeigted,1,'MULT',TS,MS);

end