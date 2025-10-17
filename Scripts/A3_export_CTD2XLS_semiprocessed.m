%Process and bin Ice2Ocean raw data

%Export  mat file with binned CTD data
%Created by  akhamilt@ualberta.ca
% modified by CP Nov. 20, 2022 to remove trimming
%May 15 2023

clearvars
cd('/Users/claire/Dropbox (Bhatia Lab)/Bhatia Lab Team Folder/Projects/Devon/Multi-year_analyses/Annual_CTD_Processed/Annual_CTD_matfiles/Andrew_RAW/');
dr= dir('*.mat');

fields = {'serialtime',...
    'temp',...
    'cond'...
    'sal',...
    'press',...
    'depth'...
    'ODO_molL',...
    'ODO_perc',...
    'turb',...
    'chla',...
    'PAR'};

for nn= 6 %1:2:length(dr)
    cd('/Users/claire/Dropbox (Bhatia Lab)/Bhatia Lab Team Folder/Projects/Devon/Multi-year_analyses/Annual_CTD_Processed/Annual_CTD_matfiles/Andrew_RAW/');
    file = dr(nn).name;  if contains(file, 'check'); continue; end;
    load(file);
    ctdraw = ctddown;
    ctd = ctddown;
    for m= 20 %1:length(ctd);
        
        %Skip over stations with no data / missing data 
        if (strcmp(ctd(m).cruiseID, 'V1908') & isnan(ctd(m).CTDserialno))    |...
            (strcmp(ctd(m).cruiseID, 'V2208') & isnan(ctd(m).CTDserialno))  
            continue
        elseif isnan(ctd(m).temp)
            continue
        end
        
        %Manually replace bad data
    	  if strcmp(ctd(m).cruiseID, "S2105") & strcmp(string(ctd(m).station), "CTDTEST1");%Cap left on Chla & ODO
                ctd(m).chla = (1:length(ctd(m).press))'*nan;
                ctd(m).ODO_molL = (1:length(ctd(m).press))'*nan;
                ctd(m).ODO_perc = (1:length(ctd(m).press))'*nan;
            elseif strcmp(ctd(m).cruiseID, "V2108") & strcmp(string(ctd(m).station), "BL15");%no measurement taken
                ctd(m).chla = (1:length(ctd(m).press))'*nan;
                ctd(m).turb = (1:length(ctd(m).press))'*nan;
          end

          %timedata = string(datetime(ctd(m).serialtime, 'convert', 'datenum'));
          % timedata(cellfun(@isempty,timedata)) = {"NaN"};

          %Trim bottom off profiles to reduce spikes
        [maxpd] = max(ctd(m).press);
        ctd(m).press(ctd(m).press>(maxpd-2)) = nan;
        ctd(m).press(ctd(m).press<0.5) = nan;   

        ctdp=ctd;
         % Low-pass filter in time (4 samples)
         %
        span = 4; 
        window = ones(span,1)/span;
            for n = 1:length(fields);
            ctdp(m).(fields{n}) = convn(ctd(m).(fields{n}),window,'valid');
            end
         %}
        % Bin ave
        ctdd=ctdp;
        ctddb = ctdd;
            if isempty(ctdd(m).press)
               disp(strcat('Error: Cannot bin ave in cast:  ', num2str(m), ' skipping to next')) 
               continue
            end
            %
            d= ctdd(m).press;
            dist = 0.25:0.25:round(max(d)); %Defines bin edges
                for n=1:length(fields);
                    ctddb(m).(fields{n}) = nan;%Fill with nan
                    x = ctdd(m).(fields{n});
                    [meanv, meand, bins]=binavgf(x,d,0.5, dist); %Number
                    ctddb(m).(fields{n}) = meanv;  
                end
            ctddb(m).press = meand;
            ctd = ctddb;
            %}
            %Derive variables
            ctd(m).SA = gsw_SA_from_SP(ctd(m).sal, ctd(m).press, ctd(m).lon, ctd(m).lat);
            ctd(m).CT = gsw_CT_from_t(ctd(m).SA, ctd(m).temp, ctd(m).press);
            ctd(m).pt = gsw_pt0_from_t(ctd(m).SA, ctd(m).temp, ctd(m).press);
            ctd(m).sigma0 = gsw_sigma0(ctd(m).SA, ctd(m).CT);
            ctd(m).rho = gsw_rho(ctd(m).SA, ctd(m).temp, ctd(m).press); % potential density
            [ctd(m).N2, ctd(m).N2press] = gsw_Nsquared(ctd(m).SA, ctd(m).temp, ctd(m).press, ctd(m).lat);
         
            O2_sol = gsw_O2sol(ctd(m).SA, ctd(m).CT, ctd(m).press, ctd(m).lon, ctd(m).lat); % umol/kg
            ctd(m).AOU =  O2_sol - ctd(m).ODO_molL*10^3*1./ctd(m).rho; % umol/kg
            ctd(m).CE = CalcCEProf(ctd(m), 'full');

    end
            newfields = fieldnames(ctd);
            for mmm = 1:length(ctd)
             for n=1:length(newfields);
                 if isempty(ctd(mmm).(newfields{n}))
                    ctd(mmm).(newfields{n}) = nan;
                 end
             end
            end

        rootfol = '/Users/claire/Dropbox (Bhatia Lab)/Bhatia Lab Team Folder/Projects/Devon/Multi-year_analyses/Annual_CTD_Processed/Annual_CTD_matfiles/Claire_PROCESSED/';
        cd(rootfol);
       filename = strcat(string(ctd(m).cruiseID), '_CTD_Data_50cm_bin.mat');
        save(filename, 'ctd')

        cd('/Users/claire/Dropbox (Bhatia Lab)/Bhatia Lab Team Folder/Projects/Devon/Multi-year_analyses/Annual_CTD_Processed/Annual_CTD_matfiles/Andrew_RAW/');
        %clearvars -except dr n m root* ctd* fields file
        %clear ctd*
end

%%
%{
figure
for n = 1:10
    plot(ctddown(n).temp, ctddown(n).press,'.-r')
    hold all
    plot(ctd(n).temp, ctd(n).press,'.-k')
end
axis ij
%}

%%

%Export csv files with binned CTD data
%Created by Claire Parrott
% Modified by akhamilt@ualberta.ca
%Mar 8 2023
%Revised to read output from RAW mat files

%Read Ice2Ocean CTD Data and output individual profiles as csv data
%Need to include units constrain sigdig, etc
%{
clearvars
cd('/Users/claire/Dropbox (Bhatia Lab)/Bhatia Lab Team Folder/Projects/Devon/Multi-year_analyses/Annual_CTD_Processed/Annual_CTD_matfiles/Andrew_RAW/');
d= dir('*.mat');
fields = {'serialtime',...
    'temp',...
    'cond'...
    'sal',...
    'press',...
    'depth'...
    'ODO_molL',...
    'ODO_perc',...
    'turb',...
    'chla',...
    'PAR'};

for nn=14%1:length(d)
    cd('/Users/claire/Dropbox (Bhatia Lab)/Bhatia Lab Team Folder/Projects/Devon/Multi-year_analyses/Annual_CTD_Processed/Annual_CTD_matfiles/Andrew_RAW/');
    file = d(nn).name;
    load(file);
    ctdraw = ctddown;
    ctd = ctddown;
    for m=1:length(ctd);    
        ctdp=ctd;
         % Low-pass filter in time (1s or 4 samples at 4Hz)
        span = 4; 
        window = ones(span,1)/span;
            for n = 1:length(fields)
            ctdp(m).(fields{n}) = convn(ctd(m).(fields{n}),window,'valid');
            end

        % Bin ave
        ctdd=ctdp;
        ctddb = ctdd;
            if isempty(ctdd(m).press)
               disp(strcat('Error: Cannot bin ave in cast:  ', num2str(m), ' skipping to next')) 
               %continue
            end
            d= ctdd(m).press;
            dist = 0:1:round(max(d));
                for n=1:length(fields);
                    ctddb(m).(fields{n}) = [];
                    x = ctdd(m).(fields{n});
                    [meanv meand bins]=binavgf(x,d,1, dist);
                    ctddb(m).(fields{n}) = meanv;  
                end
            ctddb(m).press = meand;
            ctd = ctddb;
 
            sheetdata.cruiseID= repmat(ctd(m).cruiseID, size(ctd(m).temp));
            sheetdata.station= repmat(ctd(m).station, size(ctd(m).temp));
            timedata = ones(size(ctd(m).temp))*nan;
            timedata = string(datetime(ctd(m).serialtime, 'convert', 'datenum'));
            timedata(cellfun(@isempty,timedata)) = {"NaN"};
            sheetdata.time = timedata;
            sheetdata.lat = (ctd(m).temp)*0+ctd(m).lat;
            sheetdata.lon = (ctd(m).temp)*0+ctd(m).lon;
            sheetdata.press = ctd(m).press;
            % sheetdata.depth = ctd(m).depth;
            sheetdata.temp = ctd(m).temp;
             sheetdata.cond = ctd(m).cond;
            sheetdata.sal = ctd(m).sal;
            sheetdata.ODO_molL = ctd(m).ODO_molL;
            sheetdata.ODO_perc = ctd(m).ODO_perc;
            sheetdata.turb = ctd(m).turb;
            sheetdata.chla = ctd(m).chla;
            sheetdata.PAR = ctd(m).PAR;
            %sheetdata.sigma0 = ctd(m).sigma0;
            sheetdata = struct2table(sheetdata);

% 
% figure
% for n = 2
%     plot(ctdraw(n).chla, ctdraw(n).press, 'k')
%     hold all
%     plot(ctd(n).chla, ctd(n).press, 'r')
% end
% axis ij
% 
% figure
% for n = 2
%     plot(ctdraw(n).turb, ctdraw(n).press, 'k')
%     hold all
%     plot(ctd(n).turb, ctd(n).press, 'r')
% end
% axis ij
            rootfol = '/Users/claire/Dropbox (Bhatia Lab)/Bhatia Lab Team Folder/Projects/Devon/Multi-year_analyses/Annual_CTD_Processed/Annual_CTD_csv_RAW_1mbin/';
            cruisefol = ctd(m).cruiseID;
            profiletime = ctd(m).serialtime(~isnan(ctd(m).serialtime));
            timestr = datestr(profiletime(1), 'yyyymmdd');
            if ~isdir(strcat(rootfol, cruisefol))
                mkdir(strcat(rootfol, cruisefol));
            end
            cd(strcat(rootfol, cruisefol));
            filename = strcat(string(ctd(m).cruiseID), '_', string(ctd(m).station), '_', timestr, '.csv');
            writetable(sheetdata, filename, 'FileType', 'text'); 
            clear sheetdata
       cd('/Users/claire/Dropbox (Bhatia Lab)/Bhatia Lab Team Folder/Projects/Devon/Multi-year_analyses/Annual_CTD_Processed/Annual_CTD_matfiles/Andrew_RAW/');
        clearvars -except d nn m root*  fields ctd
        end
        clear ctd
end


%}