%%%%%%%%%%%%%%%%%%%%%%%%%%
% combine nutrient data with CTD data into table called combined. Used in
% scripts that contain nutrient data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
close all

% load nutrients
cd('C:\Users\claire\Dropbox (Bhatia Lab)\Bhatia Lab Team Folder\Projects\Devon\Multi-year_analyses\Analysis\nutrients\master_sheets');
load('nutrients_multiyear_mastersheet.mat');

% filter major indicies 
idx = (nutrients.Season == 'Summer')& (nutrients.Water == 'Marine'); %&contains(nutrients.flag, 'Y');
nutrients = nutrients(idx,1:19);

% sites to use to classify riverine vs glacier sources 
sites = {'Belcher'; 'Sverdrup'; 'Jakeman'; 'Sydkap'; 'Terry'; 'Grise'; 'Fram'; 'True';'Harbour'};

open_stations = {'VIO6';'VIO7';'VIO8';'VIO32';'VIO33';
    'JS1';'JS2';'JSD';
    'JS2';'JS14';
    'JSDT_B';'OG';'S1';'JSOG1B';'SD10B'};
open_campaign = {'V1907';'V1907';'V1907';'V1907';'V1907';
    'V2007';'V2007';'V2007';
    'A2108';'V2108';
    'V2208';'V2208';'V2208';'A2208';'A2208'};

% filter by year
n2019 = nutrients(nutrients.Year == 2019,:);
n2020 = nutrients(nutrients.Year == 2020,:);
n2021 = nutrients(nutrients.Year == 2021,:);
n2022 = nutrients(nutrients.Year == 2022,:);

% loop through table, match nutrient measurement to CTD profile
cd("C:\Users\claire\Dropbox (Bhatia Lab)\Bhatia Lab Team Folder\Projects\Devon\Multi-year_analyses\Annual_CTD_Processed\Annual_CTD_matfiles\Andrew_PROCESSED");
direct = dir('*.mat');

combined = table();
errors = struct();
errors.noMatchST.stn = {'1'}; errors.noMatchST.transect = {'1'}; errors.noMatchST.transecta = {'1'}; errors.noMatchST.cruiseID = {'1'};
errors.MultipleMatches.stn = {'1'}; errors.MultipleMatches.transect = {'1'}; errors.MultipleMatches.cruiseID = {'1'};
errors.noMatch.stn = {'1'}; errors.noMatch.n = {'1'}; errors.noMatch.cruiseID = {'1'};
errors.ProfErr.stn = {'1'}; errors.ProfErr.cruiseID = {'1'}; errors.ProfErr.n = {'1'};

for n = 1:height(nutrients)
    stn = nutrients(n, :);
    cruiseID = string(stn.CruiseID);
    d = contains({direct.name}, cruiseID);
    load(direct(d).name); % load the relevant ctd data
    sid = strcmp([ctd.station], string(stn.StationID));
    if sum(sid) > 1
        %disp('uh oh, using transect to help')
        sid = strcmp([ctd.station], string(stn.StationID)) & contains([ctd.location], string(stn.Location));
        if sum(sid) > 1
            % try matching to bottom depth
            sid = strcmp([ctd.station], string(stn.StationID))&contains([ctd.location], string(stn.Location))&strcmp([ctd.bottlecast], 'Y');
            if sum(sid)>1
                sid = strcmp([ctd.station], string(stn.StationID))&contains([ctd.transect], string(stn.Transect))&strcmp([ctd.bottlecast], 'Y');
                %disp(stn.StationID); disp(cruiseID);
                if sum(sid)>1
                    errors.MultipleMatches.stn{end+1} = stn.StationID;
                    errors.MultipleMatches.cruiseID{end+1} = cruiseID;
                    errors.MultipleMatches.transect{end+1} = stn.TransectActual;
                    continue
                elseif sum(sid) == 0
                    errors.noMatchST.stn{end+1} = stn.StationID;
                    errors.noMatchST.cruiseID{end+1} = cruiseID;
                    errors.noMatchST.transect{end+1} = stn.TransectActual;
                end
            elseif sum(sid) == 0
                    %disp('uh oh, nothing matches still!')
                    %disp(stn.StationID); disp(cruiseID);
                    errors.noMatchST.stn{end+1} = stn.StationID;
                    errors.noMatchST.cruiseID{end+1} = cruiseID;
                    errors.noMatchST.transect{end+1} = stn.TransectActual;
                    continue
            end
        elseif sid == 0
            % the location match didn't work.. bottlecast only
            if stn.StationID == 'JSOG1B'
                sid = strcmp([ctd.station], string(stn.StationID))&strcmp([ctd.bottlecast], 'N');
            elseif stn.StationID == 'JSDT_B'
                sid = strcmp([ctd.station], 'JSDT') & strcmp([ctd.bottlecast], 'Y');
                sid(42) = 0;
            else
                sid = strcmp([ctd.station], string(stn.StationID)) & strcmp([ctd.bottlecast], 'Y');
            end
            if sum(sid)>1
                %disp('bigger problems');
                %disp(stn.StationID); disp(cruiseID);
                errors.MultipleMatches.stn{end+1} = stn.StationID;
                errors.MultipleMatches.cruiseID{end+1} = cruiseID;
                errors.MultipleMatches.transect{end+1} = stn.TransectActual;
                continue
            elseif sum(sid) == 0
                %disp('uh oh, nothing matches still!')
                %disp(stn.StationID); disp(cruiseID);
                errors.noMatchST.stn{end+1} = stn.StationID;
                errors.noMatchST.cruiseID{end+1} = cruiseID;
                errors.noMatchST.transect{end+1} = stn.TransectActual;
                continue
            end
        end
    elseif sid == 0
            if stn.StationID == 'JSDT_B'
                sid = strcmp([ctd.station], 'JSDT') & strcmp([ctd.bottlecast], 'Y');
                sid(42) = 0;
            end
    end
  %
    if sum(sid) == 0
        if any(nutrients(n,:).StationID == {'SSW', 'DSW'})
            continue
        else
            %disp('uh oh, nothing matches?')
            %disp(stn.StationID); disp(cruiseID);
            errors.noMatch.stn{end+1} = stn.StationID;
            errors.noMatch.n{end+1} = n;
            errors.noMatch.cruiseID{end+1} = cruiseID;
            disp(n);
        continue
        end
    end

    %% label location of open stations
    if any(stn.StationID == open_stations)% && any(stn.CruiseID == open_campaign);
         stn.Location = 'Jones Sound';
    end

    profile = ctd(sid);
    if isnan(profile.temp); disp(profile.station); continue; end
    % distance from coast
    [dist, ~] = find_distance(profile);
    %% add physical metrics
    [dz] = Getdz(profile.station, profile.location, profile.cruiseID);
    %{
    if isnan(dz) && strcmp(profile.transect, 'Jakeman Glacier - practice')
        [dz] = GetUpwellingIndex(profile.station, 'Jakeman practice', profile.cruiseID);
    end
    %}
    if isnan(dz) && (any((stn.Location ~= 'Jones Sound')))
        dz = 0;
    end
    %% calculate full column freshwater metrics 
    [maxN2, ~]  = CalcMaxN2(profile);

    %{
    if stn.Depth_m <= 30
        [FWT, ~]    = CalcFWT(profile, 33., 100);
    elseif stn.Depth_m > 30 | stn.Depth_m <= 250
        [FWTs, ~]    = CalcFWT(profile, 33., 10);
        [FWTi, ~]    = CalcFWT(profile, 33., 250);
        FWT = FWTi - FWTs;
    else
        [FWT, ~]    = CalcFWT(profile, 34.8, []);
    end
    %}
    [FWT, ~]    = CalcFWT(profile, 34.8, 100);
    [CE100, ~,~,~] = CalcCE(profile, 100);
    [CE30, ~,~,~] = CalcCE(profile, 30);
    [CErh27, ~,~,~] = CalcCE(profile, [],'rho27');
    [CErh26, ~,~,~] = CalcCE(profile, [],'rho26');
    [CEd,~,~,~] = CalcCE(profile, stn.Depth_m);
    %% Added Chem metrics
    O2_sol = gsw_O2sol(profile.SA, profile.CT,profile.press, profile.lon, profile.lat); % umol/kg
    AOU = O2_sol - profile.ODO_molL*10^3*1./profile.rho; % umol/kg

    %%
    [~,didx] = min(abs(profile.depth - stn.Depth_m)); % get the profile matching the station and the measurement at that depth
    if didx == 1
        %disp('issue with profile')
        %disp(stn.StationID); disp(cruiseID);
        errors.ProfErr.stn{end+1} = stn.StationID;
        errors.ProfErr.n{end+1} = n;
        errors.ProfErr.cruiseID{end+1} = cruiseID;
        continue
    end
    T = stn; 

    %
    % add site classifications
    ri_idx = logical(contains(string([T.Transect]), {sites{6:9}}));
    fj_idx = logical(contains(string([T.Transect]), {sites{4:7}}));
    op_idx = logical(contains(string([T.Location]), {'OG';'Jones Sound'}));
    FWsource = repmat({'Glacier'}, [height(T), 1]);
    FWsource(ri_idx) = {'River'};
    FWsource(op_idx) = {'Open'};
    coast = repmat({'Open'}, [height(T), 1]);
    coast(fj_idx) = {'Fjord'};
    coast(op_idx) = {'JS'};
    for mm = 1:height(T)
    comb_char{mm} = char(FWsource{mm} + "-" + coast{mm});
    comb_chars{mm} = char(FWsource{mm}(1) + "-" + coast{mm}(1));
    end
    
    %% glacier melt 
    if strcmp(FWsource, 'Glacier')
        idg = find(contains(string([T.Transect]), {sites{1:4}}));
        if isempty(idg);   glmelt = nan;        else;
        gl_site = sites{idg}; % 
        load(strcat('C:\Users\claire\Dropbox (Bhatia Lab)\Bhatia Lab Team Folder\Projects\Devon\Multi-year_analyses\Mass_Balance\Bhatia\Output_matfiles\',gl_site, '-DailyRACMOOutput.mat'));
        % get discharge to day of sampling
        tmask = glacier_struct.time >= datetime(T.Year, 1,1) & glacier_struct.time <= T.Date_yyyy_mm_dd_;
        tmask1 = glacier_struct.time == T.Date_yyyy_mm_dd_;
        glmelt = sum(glacier_struct.runoff(tmask), 'omitnan'); % kg/m3 ?
        glmelt1 = sum(glacier_struct.runoff(tmask1), 'omitnan'); % kg/m3 ?
        end
    else
        glmelt = nan;
        glmelt1 = nan;
    end
    T.ctd_name = profile.station;
    %
    T.FWSource = categorical(FWsource); T.Coastline = categorical(coast); T.combinedchar = categorical(comb_char); T.combinedcharshort = categorical(comb_chars);
    T.dist = round(dist, 1);
    T.FWT = FWT; T.CE =  CE100; T.CE100 = CE100 - CE30; T.CE30 = CE30; T.CE = CEd; T.CErh = CErh26; T.CErh = CErh27; T.maxN2 = maxN2;
    T.pt = profile.pt(didx); T.CT = profile.CT(didx); T.sa = profile.SA(didx); T.rho = profile.rho(didx); 
    T.DO_Lmol = profile.ODO_molL(didx); T.DO_perc = profile.ODO_perc(didx); 
    T.chla = profile.chla(didx); T.turb = profile.turb(didx); T.sigma = profile.sigma0(didx);
    T.dz = dz; T.AOU = AOU(didx); T.glmelt = glmelt; T.glmelt1 = glmelt1;
    %% get depth average for pertinent variables
    d_30 = find(profile.depth >= 30.00, 1); 
    T.turb_c = mean(profile.turb(1:d_30)); T.ODO_c = mean(profile.ODO_molL(1:d_30)); T.chla_c = mean(profile.chla(1:d_30));  T.sa_c = mean(profile.SA(1:d_30));  T.pt_c = mean(profile.pt(1:d_30));
    d_100 = find(profile.depth >= 100.00, 1); if isempty(d_100); d_100 = length(profile); end;
    T.turb_d = mean(profile.turb(d_30:d_100)); T.ODO_d = mean(profile.ODO_molL(d_30:d_100)); T.chla_d = mean(profile.chla(d_30:d_100));  T.sa_d = mean(profile.SA(d_30:d_100));  T.pt_c = mean(profile.pt(1:d_30));
    T.loopnumber = n;
    if any([isnan(T.turb_c), isnan(T.ODO_c), isnan(T.chla_c)])
        mask = ~isnan(profile.depth); d_30 = find(profile.depth(mask) >= 30.00, 1);
        p_turb = profile.turb(mask);    T.turb_c = mean(p_turb(1:d_30));
        p_ODO = profile.ODO_molL(mask); T.ODO_c = mean(p_ODO(1:d_30));
        p_chla = profile.chla(mask);    T.chla_c = mean(p_chla(1:d_30));
    end

    %{
    %% something with PAR. Not certain about these units but assumuing this is percent of surface light available?
    if any(strcmp(fields(profile), 'PAR')); if ~isempty(profile.PAR);
        id_pd = find(profile.PAR < 0.01, 1); T.PARdepth = profile.depth(id_pd);
    else T.PARdepth = nan;
    end; end
    %}  
    %% 
    if didx == length(profile.depth)
        T.N2 = profile.N2(didx-1);
    else
        T.N2 = profile.N2(didx);
    end
    combined(n,:) = T;
    clear T didx sid AOU CE* FWT 
end

combined = combined(~(combined.Year == 0),:); % save only non-empty entries

cd('C:\Users\claire\Dropbox (Bhatia Lab)\Bhatia Lab Team Folder\Projects\Devon\Multi-year_analyses\Analysis\nutrients\master_sheets\Combined_CTD')
save('combined_nuts_ctd_all_noflag_dz5.mat', 'combined');


