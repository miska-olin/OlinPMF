function [ out ] = readHrSpectra( filename, path )
%READHRSPECTRA Use this function to import HR spectra preprocessed in tofTools


    % Check input arguments
    if nargin < 1
        error('Too few input arguments')
    elseif nargin == 1
        filepath = filename;
    elseif nargin == 2
        filepath = [path filename];
    elseif nargin > 2
        error('Too many input arguments')
    end
    
    % Load UMR stick data vector
    % data = h5read(filepath,'/Sticks/signal');
    % data = data';

    % Load time vector
    time = h5read(filepath,'/tofData/time');

    % Load peaklist
    peaklist = h5read(filepath, '/GUIdata/peakListCompounds');
    peakdatagauss = h5read(filepath, '/Peaks/PeakDataGauss');
    peaklistselected = h5read(filepath, '/GUIdata/peakListSelected');

    if ~all(peaklistselected)
        warning('HR peak time series can only be exported with selecting all peaks from the peak list. Peak fitting must also be performed after selecting them.')
    
        peaks = struct;
    else
        
        peaks.compound = {};
        peaks.HR = single(nan(length(time),length(peaklist)));

        disp("peakdatagauss size "+size(peakdatagauss,1)+"x"+size(peakdatagauss,2) + " peaklist length "+length(peaklist) + " time length " + length(time))
        disp("   " + size(peakdatagauss,2) + "/"+length(peaklist) + " = " + size(peakdatagauss,2)/length(peaklist))
        for i=1:length(peaklist)
            if mod(i,10)==1
                hw=waitbar(i/length(peaklist));
            end

            peaks.compound{i} = deblank(peaklist{i});
            if isnan(str2double(peaks.compound{i}))
                if exist('tof_exact_mass.m','file')==0
                    error('tofTools package is needed to initialize variables from a tofTools h5 file')
                end
                peaks.mass(i) = tof_exact_mass(peaks.compound{i});
            else
                peaks.mass(i) = str2double(peaks.compound{i});
            end
            peaks.UMR(i) = round(peaks.mass(i));

            I=find(peakdatagauss(8,:)==i);
            % disp("numel I= " + numel(I))
            try
                timeIndexes = peakdatagauss(1,I);
                peaks.HR(timeIndexes,i)=peakdatagauss(5,I);
            catch ME
                disp("i="+i+ " length(I)="+length(I) + " compound: " + peaks.compound{i})
                rethrow(ME)
            end
        end

        close(hw)

    end
    
        
    peaks.elements = array2table(h5read(filepath, '/Peaks/PeakElements')');
    elementnames = h5readatt(filepath, '/Peaks/PeakElements','Elements');
    elementnames = split(elementnames);
    for i=1:length(elementnames)-1
        peaks.elements.Properties.VariableNames{i} = elementnames{i};
    end
    
    % % Load acquisition log
    % acq_log = h5read(filepath,'/Logs/AcquisitionLog');
    % log = struct;
    % log.note = {};
    % for i=1:length(acq_log.Time)
    %     log.timestamp(i) = datetime(acq_log.Time(i),'ConvertFrom','datenum');
    %     log.note{i} = acq_log.Message{i};
    % end
    
    out = struct;
    % out.data = data;
    out.time = datetime(time','ConvertFrom','datenum');
    out.peaks = peaks;
    % out.log = log;
end

