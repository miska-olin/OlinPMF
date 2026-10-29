
hf=figure(1);
hf.Position=[1351 19 560 784];

subplot(211)
yyaxis left
plot(p_vec,QQexp_robust_min_vec,'x-',p_vec,ones(size(p_vec)),'-k')
ylabel('Q/Qexp robust')

yyaxis right
plot(p_vec,QQexp_true_min_vec,'x-',p_vec,ones(size(p_vec)),'-k')
xlabel('p')
ylabel('Q/Qexp true')


title('PMF computation running.\newlineClose this window to abort computation while keeping the current result.',FontWeight='normal',FontSize=8)



subplot(212)
plot(p_vec,100*unExplainedVariation_vec,'ok-',p_vec,zeros(size(p_vec)),'-k')
ylabel('Unexplained fraction (%)')
xlabel('p')

drawnow