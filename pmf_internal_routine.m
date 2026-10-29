function [TS,MS,dnorm] = pmf_internal_routine(X,TS0,MS0,alg,S,nDownWeighted)
    
max_i=100;
    for i=1:max_i

        if strcmp('MULT',alg)
            temp = TS0'*X;
            MS = max(0, MS0.*(temp./((TS0'*TS0)*MS0 + eps(temp))));

            temp = X*MS';
            TS = max(0, TS0.*(temp./(TS0*(MS*MS') + eps(temp))));
        else

            MS = max(0, TS0\X);
            TS = max(0, X/MS);
        end
        
        EPerS=(X-TS*MS)./S;
        EPerS(abs(EPerS)>4)=4;
    
        dnorm=(sum(sum((EPerS).^2))/(numel(X)-nDownWeighted-numel(TS)-numel(MS)));
             
    
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

