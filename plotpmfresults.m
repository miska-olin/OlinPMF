if useGPU
    TS_final=gather(TS_final);
    MS_final=gather(MS_final);
    S=gather(S);
    X_meas=gather(X_meas);
    X_meas_nonConstrained=gather(X_meas_nonConstrained);
    mzs=gather(mzs);
end


if nConstrainedFactors==0
    E_final=X_meas-TS_final*MS_final;
else
    if involveConstrainedFactorsInQCalc==1
        MS_final_con_full=[MS_final; constrainedMS];
        TS_final_con_full=[TS_final constrainedTS];
        E_final=X_meas_nonConstrained-TS_final_con_full*MS_final_con_full;

    else
        E_final=X_meas-TS_final*MS_final;
    end
end



hf=figure(p);
hf.Position=[97 1 1236 994];

cmap=colormap(hsv(p));
clear ax

for i=1:p


    ax(i)=subplot(p,2,i*2-1);
    if forceUMR
        bar(round(mzs),MS_final(i,:),'BarWidth',1,'EdgeColor','k','facecolor',cmap(i,:))
    else
        if ~strcmp(pmfType,'HR')
            bar(mzs+binWidth/2,MS_final(i,:),'BarWidth',1,'EdgeColor','k','facecolor',cmap(i,:))
        else % HR
            hbar=bar(mzs,MS_final(i,:),'BarWidth',10,'EdgeColor','flat','facecolor',cmap(i,:));

            if printMassTexts
                if printOnlyHighestMassTexts==0
                    for iMass=1:length(mzs)
                        text(mzs(iMass)+2,MS_final(i,iMass),compoundNames(iMass),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',7,'Rotation',0)
                    end
                else
                    [~,I]=sort(MS_final(i,:));
                    for iMass=I(end-printOnlyHighestMassTexts:end)
                        text(mzs(iMass)+2,MS_final(i,iMass),compoundNames(iMass),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',7,'Rotation',0)
                    end
        
                end
            end
        
            if addMassTextDataTips
                hbar.DataTipTemplate.DataTipRows(3)=dataTipTextRow('Mass Text',compoundNames);
            end
        end
    end
    
    ylabel(['F' int2str(i)])

end

xlabel('m/z (Th)')
linkaxes(ax,'x')

iCutPoint=ceil(p/2)*2;




if plotUnexpAsTime

    axTS(1)=subplot(p,2,2:2:iCutPoint);


    yyaxis left
        denominator=zeros(size(TS_final,1),1);
        for jP=1:length(denominator)
            for iP=1:size(S,2)
                temp = 0;
                for hP=1:(p+nConstrainedFactors*involveConstrainedFactorsInQCalc)
                    if nConstrainedFactors==0
                        temp=temp+MS_final(hP,iP)*TS_final(jP,hP);
                    else
                        temp=temp+MS_final_con_full(hP,iP)*TS_final_con_full(jP,hP);
                    end
                end
        
                denominator(jP)=denominator(jP)+(temp+abs(E_final(jP,iP)))/S(jP,iP);
            end
        end
        numerator=sum(abs(E_final)./S,2);
        plot(datetime(t,'convertfrom','datenum'),100*numerator./denominator);
        xlabel('Time')
        ylabel('Unexplained variation (%)')
        
        unExplainedVariation=sum(numerator)/sum(denominator);
        title({['Q/Qexp = ' num2str(QQexp_final)  ' (robust: ' num2str(QQexp_robust_min) ', true: '  num2str(QQexp_true_min) ')'],['Unexplained variation = ' num2str( unExplainedVariation*100,3 ) '%']})

     yyaxis right
        plot(datetime(t,'ConvertFrom','datenum'),(sum(min((E_final./S).^2,inf),2))/size(MS_final,2))
        ylabel('Q/Qexp(t)')


else

    subplot(p,2,2:2:iCutPoint)

    denominator=zeros(size(MS_final,2),1);
    for jP=1:length(denominator)
        for iP=1:size(S,1)
            temp = 0;
            for hP=1:(p+nConstrainedFactors*involveConstrainedFactorsInQCalc)
                if nConstrainedFactors==0
                    temp=temp+TS_final(iP,hP)*MS_final(hP,jP);
                else
                    temp=temp+TS_final_con_full(iP,hP)*MS_final_con_full(hP,jP);
                end
            end
    
            denominator(jP)=denominator(jP)+(temp+abs(E_final(iP,jP)))/S(iP,jP);
        end
    end
    numerator=sum(abs(E_final)./S,1);
    bar(mzs+binWidth/2,100*numerator./denominator','BarWidth',1,'EdgeColor','k')
    xlabel('m/z (Th)')
    ylabel('Unexplained variation (%)')
    
    unExplainedVariation=sum(numerator)/sum(denominator);
    title({['Q/Qexp = ' num2str(QQexp_final)  ' (robust: ' num2str(QQexp_robust_min) ', true: '  num2str(QQexp_true_min) ')'],['Unexplained variation = ' num2str( unExplainedVariation*100,3 ) '%']})




end


if plotUnexpAsTime
    axTS(2)=subplot(p,2,(iCutPoint+2):2:(2*p));
else
    subplot(p,2,(iCutPoint+2):2:(2*p))
end

legendNames={''};


for i=1:p


    
    plot(datetime(t,'convertfrom','datenum'),TS_final(:,i),'color',cmap(i,:));

    hold on

    legendNames{i}=['F' int2str(i)];
end

hold off

if plotUnexpAsTime
    linkaxes(axTS,'x')
end

legend(legendNames)
ylabel('Factor Time Series')


drawnow

zoom on


if useGPU
    TS_final=gpuArray(TS_final);
    MS_final=gpuArray(MS_final);
    S=gpuArray(S);
    X_meas=gpuArray(X_meas);
    X_meas_nonConstrained=gpuArray(X_meas_nonConstrained);
    mzs=gpuArray(mzs);
end