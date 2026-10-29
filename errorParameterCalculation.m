% USE THIS AFTER RUNNING A QUICK PMF USING THE CORRECT TYPE (UMR, BIN, HR)
% FIRST SO THAT THE VARIABLES ARE LOADED CORRECTLY


nDecile=10; %advanced parameter, default=10

sigma=nan(size(X_meas,2),nDecile);
sigmaError=sigma;
meanCps=sigma;


for i=1:length(mzs)
    
    mi=double(X_meas(:,i));
    miSorted=sort(mi);

    nSorted=sum(~isnan(miSorted));
    
    decileLimits=(miSorted(1:ceil(nSorted/10):end));
    mi(mi>decileLimits(10))=nan;
    miSorted2=sort(mi);
    nSorted2=sum(~isnan(miSorted2));
    
    mimedian=movmedian(mi,5);
    mimedianSorted=sort(mimedian);
    nmedianSorted=sum(~isnan(mimedianSorted));
    decileLimits2=mimedianSorted(1:ceil(nmedianSorted/nDecile):nmedianSorted);
    
    diff_i=mi-mimedian;
    diff_i(abs(diff_i)==0)=nan;

    for iDecile=1:nDecile

        if iDecile==nDecile
            I=mimedian>decileLimits2(iDecile);
        else
            I=(mimedian<decileLimits2(iDecile+1)&mimedian>decileLimits2(iDecile));
        end

        
        meanCps(i,iDecile)=(mean(sqrt(mi(I&mi>=0)),'omitnan'))^2;


        h=histogram(diff_i(I),'BinMethod','scott');
        
        try
            [fitoutput,gof,output]=fit((h.BinEdges(1:end-1)+h.BinEdges(2:end))'/2,h.BinCounts','gauss1');
            fitoutputError=nlparci([fitoutput.a1 fitoutput.b1 fitoutput.c1],output.residuals,'jacobian',output.Jacobian,'alpha',0.3174);

            sigma(i,iDecile)=fitoutput.c1/sqrt(2);
            sigmaError(i,iDecile)=mean(abs(fitoutputError(3,:)-fitoutput.c1))/sqrt(2);

            
            if abs(fitoutput.b1)>fitoutput.c1/sqrt(2)/2||abs(mean(h.BinEdges(1:2)))<1.0*fitoutput.c1/sqrt(2)||abs(mean(h.BinEdges((end-1):end)))<1.0*fitoutput.c1/sqrt(2)
                sigma(i,iDecile)=nan;
                sigmaError(i,iDecile)=nan;
            end



        catch
            sigma(i,iDecile)=nan;
            sigmaError(i,iDecile)=nan;
        end
        

    end
end

close all

curveFitter(meanCps,sigma,[],sigma./(sigmaError.^2))
disp(' ')
disp('Fit a custom equation using fitting parameters you choose, after removing clear outliers.')
disp('A good idea is to use the following function as a custom equation (which is the same form as runpmf.m assumes:')
disp('sqrt(errorParameter1^2*x+errorParameter2^2)')
disp(' ')
disp('Note that both errorParameters should be positive (errorParameter2 could be very close to 0), and that they depend on the averaging time of the data.')
disp('After obtaining the errorParameters for the specific PMF type (UMR, bin, HR), insert them into runpmf.m settings.')