% Loads relevant data from Metadata log sheet and creates .mat structures of transects indicated by the log.
%%% Input is .rsk file, only one file can be processed at a time. Script is
%%% run as many times as there are log sheets. Uses processrsk.m to process
%%% profiles. 

%Run 'A1_read_metadata_csv2mat_dynamic.m' before this code
%Plot the output from this code using 'A2b_plot_raw_CTD_data.m'

%Created sometime in 2022
% Written by Andrew Hamilton (akhamilt@ualberta.ca) & Claire Parrott
%(cparrott@eoas.ubc.ca)
%Last modified: Dec 2 2022

clearvars
close all

%Load metadata log sheets from mat format
cd('/Users/claire/Dropbox (Bhatia Lab)/Bhatia Lab Team Folder/Projects/Devon/Multi-year_analyses/Annual_Logs_Standardized/Annual_CTD_Logs_Standardized_mat/');
d = dir('*.mat');
uberfields = {'cruiseID', 'location', 'station', 'date', 'time', 'lon', 'lat', 'bottom_depth', 'bottlecast', 'transect', 'metadata', 'bottle_depths',...
'CTDserialno', 'serialtime', 'press', 'temp', 'cond', 'sal', 'turb', 'chla', 'ODO_molL', 'ODO_perc', 'depth', 'PAR'};

for n = 12 %[1 3 12 14 15 16] %1:length(d)
    logfilename = d(n).name
    load(logfilename);
    %Find folder with relevant CTD RSK files based on cruise name
    %Get cruise name
    refcruiseID = logfilename(1:5)
    %Change name of cruise specific metadata table to generic 'log'
    log = metadata;
    %CD to root folder
    ctdroot = '/Users/claire/Dropbox (Bhatia Lab)/Bhatia Lab Team Folder/Projects/Devon/Multi-year_analyses/Annual_CTD_Data_Raw/';
    cd([ctdroot refcruiseID '_CTD_Data/'])
    %Run through each line in log and read associated rsk file
    rskfile = []; nn = 1; 
    %Create struct
    ctd = struct();
    for kk = 1:length(uberfields)
        ctd.(uberfields{kk}) = [];
        ctddown = ctd;
        ctdup = ctd;
    end

    emptyprofs = zeros(1, length(log.CruiseID));
    for nn = 42 %1:length(log.CruiseID)
         %Read CTD rsk file
        rskfile = char(log.OriginalRSKfilename(nn));
        profileno = log.RSKFullProfileno(nn);
        if isnan(profileno)
            emptyprofs(nn) = 1;
        end
        %if profileno>=1
        %Populate metadata
        ctd(nn).cruiseID = log.CruiseID(nn);
        ctd(nn).location = log.Location(nn);
        ctd(nn).station = log.StationID(nn);
        ctd(nn).date = log.Dateyyyymmdd(nn);%
        ctd(nn).time = log.StartTimeHHMM(nn);%
        ctd(nn).lon = log.LonddddddE(nn);
        ctd(nn).lat = log.LatddddddN(nn);
        ctd(nn).bottom_depth = log.BottomDepthm(nn);
        ctd(nn).bottlecast = log.BottleCast(nn);
        ctd(nn).transect = log.Transect(nn);
        ctd(nn).metadata = log(nn,:);

        %Read in bottle depths
        varnames = log.Properties.VariableNames;
        bb=1;
        ctd(nn).bottle_depths = nan;
        for xx = 1:length(log.Properties.VariableNames)
            if     contains(varnames{xx}, 'Bottle_') && contains(varnames{xx}, 'Depth_m')
             ctd(nn).bottle_depths(bb) = double(string(log.(varnames{xx})(nn)));
             bb = bb+1;
            end
        end
        %Remove extra nans
        ctd(nn).bottle_depths = ctd(nn).bottle_depths(~isnan(ctd(nn).bottle_depths));
       ctd(nn).bottle_depths(isempty(ctd(nn).bottle_depths)) = nan;

  
        %Open new rsk file
        %Read data from rsk file
        for m = 1:2 %This is to read downcast then upcast
            %FIll in field names for metadata entries without an associated CTD profile
            if isnan(profileno);                        
                fields = uberfields;
                for jj = 1:length(fields)
                    if isempty(ctd(nn).(fields{jj}))
                        ctd(nn).(fields{jj}) = nan;
                    end
                end
            else
            rsk = RSKopen(rskfile);
            rskt = RSKreaddata(rsk);
            
            %Correction for S2212 cruise where CTD sampling interval too
            %low (12s) so no casts detected in RUSKIN
            if strcmp(rskfile, '206288_20221222_1948 December 2022.rsk');
            rsk2 = RSKtimeseries2profiles(rskt);
            rskt=[];
            rskt = rsk2;
            end
            
            %Correction for V2208 where errors in simulacast rsk file so no
            %casts detected in Ruskin
            if strcmp(rskfile, '201462_20220804_2258_DOtest.rsk');
            rsk2 = RSKtimeseries2profiles(rskt);
            rskt=[];
            rskt = rsk2;
            end

            % read the downcast from profile     
                if m==1
                 rsk = RSKreadprofiles(rskt,'profile',profileno,'direction', 'down');
                elseif m==2
                 rsk = RSKreadprofiles(rskt,'profile',profileno,'direction', 'up');
                end
            %Convert channels names to cell vector and trim trailing whitespace
            for mm = 1:length(rsk.channels);
            rsk.channels(mm).longName = strtrim(rsk.channels(mm).longName);
            end
            %Convert data format to my struct format
            ctd(nn).CTDserialno = num2str(rsk.instruments.serialID);
            ctd(nn).serialtime = rsk.data.tstamp; %should convert all times to UTC for consistency
            ycol = getchannelindex(rsk, 'Sea Pressure');
            ctd(nn).press      = rsk.data.values(:,ycol);
            ycol = getchannelindex(rsk, 'Temperature');
            ctd(nn).temp       = rsk.data.values(:,ycol);
            ycol = getchannelindex(rsk, 'Conductivity');
            ctd(nn).cond      = rsk.data.values(:,ycol);
            ycol = getchannelindex(rsk, 'Salinity');
            ctd(nn).sal      = rsk.data.values(:,ycol);
            ycol = getchannelindex(rsk, 'Turbidity');
            ctd(nn).turb         = rsk.data.values(:,ycol);
            ycol = getchannelindex(rsk, 'Chlorophyll a');
            ctd(nn).chla     = rsk.data.values(:,ycol);
            ycol = getchannelindex(rsk, 'Dissolved O21');
            ctd(nn).ODO_molL     = rsk.data.values(:,ycol);
            ycol = getchannelindex(rsk, 'Dissolved O22');
            ctd(nn).ODO_perc     = rsk.data.values(:,ycol);
            ycol = getchannelindex(rsk, 'Depth');
            ctd(nn).depth     = rsk.data.values(:,ycol);
            if rsk.instruments.serialID== 201462 ||  (rsk.instruments.serialID == 206288 && ctd(nn).serialtime(1) > datenum('2022-03', 'yyyy-mm'));%Only Maya's CTD (201462) had PAR initiall. PAR installed on Luke's (206288) in June 2022.
            ycol = getchannelindex(rsk, 'PAR');
            ctd(nn).PAR               = rsk.data.values(:,ycol);
            else 
                ctd(nn).PAR = ((1:length(ctd(nn).temp))*nan)';%Create PAR field with NaNs if no PAR sensor
            end
    
            clear rsk
            end

      if m==1;
            ctddown(nn) = ctd(nn);
        elseif m==2
            ctdup(nn) = ctd(nn);
      end
        end
    	nn = nn+1;
        end

    %ctddown = ctddown(~logical(emptyprofs));
    %ctdup = ctdup(~logical(emptyprofs));
    %{
    cd('/Users/claire/Dropbox (Bhatia Lab)/Bhatia Lab Team Folder/Projects/Devon/Multi-year_analyses/Annual_CTD_Processed/Annual_CTD_matfiles/Andrew_RAW/')
    file2save = string([refcruiseID,'_CTD_Data_RAW_downup-check.mat']);
    save(file2save, 'ctddown', 'ctdup')
    clearvars -except n nn d uberfields
    cd('/Users/claire/Dropbox (Bhatia Lab)/Bhatia Lab Team Folder/Projects/Devon/Multi-year_analyses/Annual_Logs_Standardized/Annual_CTD_Logs_Standardized_mat/')
    %}
end


