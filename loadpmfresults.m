function loadpmfresults(nFactors,type)

switch type
    case 'bin'
        
        evalin('base',['p = ' int2str(nFactors) ';'])
        evalin('base',['MS_final = binPMFRes' int2str(nFactors) 'Factors.MS;'])
        evalin('base',['TS_final = binPMFRes' int2str(nFactors) 'Factors.TS;'])
        evalin('base','X_meas = X_meas_bin;')
        evalin('base','S = S_bin;')
        evalin('base','mzs = mzs_bin;')
        evalin('base','pmfType = ''bin'';')
  
        evalin('base','if useGPU,TS_final=gather(TS_final);,end')
        evalin('base','if useGPU,MS_final=gather(MS_final);,end')
        evalin('base','if useGPU,X_meas=gather(X_meas);,end')


    case {'umr','UMR'}

        evalin('base',['p = ' int2str(nFactors) ';'])
        evalin('base',['MS_final = UMRPMFRes' int2str(nFactors) 'Factors.MS;'])
        evalin('base',['TS_final = UMRPMFRes' int2str(nFactors) 'Factors.TS;'])
        evalin('base','X_meas = X_meas_UMR;')
        evalin('base','S = S_UMR;')
        evalin('base','mzs = mzs_UMR;')
        evalin('base','pmfType = ''UMR'';')


        evalin('base','if useGPU,TS_final=gather(TS_final);,end')
        evalin('base','if useGPU,MS_final=gather(MS_final);,end')
        evalin('base','if useGPU,X_meas=gather(X_meas);,end')

            

    case {'hr','HR'}

        evalin('base',['p = ' int2str(nFactors) ';'])
        evalin('base',['MS_final = HRPMFRes' int2str(nFactors) 'Factors.MS;'])
        evalin('base',['TS_final = HRPMFRes' int2str(nFactors) 'Factors.TS;'])
        evalin('base','X_meas = X_meas_HR;')
        evalin('base','S = S_HR;')
        evalin('base','mzs = mzs_HR;')
        evalin('base','pmfType = ''HR'';')


        evalin('base','if useGPU,TS_final=gather(TS_final);,end')
        evalin('base','if useGPU,MS_final=gather(MS_final);,end')
        evalin('base','if useGPU,X_meas=gather(X_meas);,end')


end


