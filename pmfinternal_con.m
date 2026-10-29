function [TS_best,MS_best,dnorm_best] = pmfinternal_con(X,p,S,nDownWeighted,nReps,pmfalg,TS0,MS0,X_full,nConstrainedFactors,involveConstrainedFactorsInQCalc,constrainedTS,constrainedMS)

iterationloop = @iterationLoop;

temp0=inf;
for i = 1:nReps
    temp1 = iterationloop(i);
    
    if temp1{1}<temp0
        temp0=temp1{1};
        bestResult=temp1;
    end
end

dnorm_best = bestResult{1};
TS_best = bestResult{3};
MS_best = bestResult{4};

    function out = iterationLoop(i)
    
        temp = cell(4,1);
    
        if ~isempty(TS0) && i==1  
            temp{3} = TS0;
        else
            temp{3} = rand(size(X,1),p); 
        end

        if ~isempty(MS0) && i==1 
            temp{4} = MS0;
        else
            temp{4} = rand(p,size(X,2)); 
        end
    
        [temp{3},temp{4},temp{1}] = pmf_internal_routine_con(X,temp{3},temp{4},pmfalg,S,nDownWeighted,X_full,nConstrainedFactors,involveConstrainedFactorsInQCalc,constrainedTS,constrainedMS);
        temp{2} = i;
    
        out = temp;
    end

end  
