function [TS,MS,dnorm] = pmf_internal_routine_con(X,TS0,MS0,alg,S,nDownWeighted,X_full,nConstrainedFactors,involveConstrainedFactorsInQCalc,constrainedTS,constrainedMS)
    
max_i=100;
    for i=1:max_i

        if strcmp('mult',alg)
            temp = TS0'*X;
            MS = max(0, MS0.*(temp./((TS0'*TS0)*MS0 + eps(temp))));

            temp = X*MS';
            TS = max(0, TS0.*(temp./(TS0*(MS*MS') + eps(temp))));
        else

            MS = max(0, TS0\X);
            TS = max(0, X/MS);
        end

        if involveConstrainedFactorsInQCalc == 1
            d=X_full-[TS constrainedTS]*[MS; constrainedMS];
        else
            d = X - TS*MS;
        end
        
        EPerS=d./S;
        EPerS(abs(EPerS)>4)=4;
    
        if involveConstrainedFactorsInQCalc == 1
            dnorm=sum(sum((EPerS).^2))/(numel(X)-nDownWeighted-(size(TS,2)+nConstrainedFactors)*(size(X_full,1)+size(X_full,2)));
        else
            dnorm=sum(sum((EPerS).^2))/(numel(X)-nDownWeighted-(size(TS,2))*(size(X,1)+size(X,2)));
        end
             
        dTS = max(max(abs(TS-TS0) / (sqrt(eps)+max(max(abs(TS0))))));
        dMS = max(max(abs(MS-MS0) / (sqrt(eps)+max(max(abs(MS0))))));

        if i>1
            if max(dTS,dMS) <= 1e-4
                break;

            elseif (dnorm0-dnorm) <= 1e-4*max(1,dnorm0)
                break;

            elseif i==max_i
                break;

            end
        end
    
        dnorm0 = dnorm;
        TS0 = TS;
        MS0 = MS;
    end

end

