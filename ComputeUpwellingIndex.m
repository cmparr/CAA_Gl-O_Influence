%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute upwelling index using CTD data. Index and its constituents
% (isopycnal uplift and anomalies) are stored and saved in a table.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear

% Set ctd directory and load combined dataset
ctd_root = 'C:\Users\claire\Dropbox (Bhatia Lab)\Bhatia Lab Team Folder\Projects\Devon\Multi-year_analyses\Annual_CTD_Processed\Annual_CTD_matfiles\Claire_PROCESSED\';
dts = dir([ctd_root '*.mat']);
load('C:\Users\claire\Dropbox (Bhatia Lab)\Bhatia Lab Team Folder\Projects\Devon\Multi-year_analyses\Analysis\nutrients\master_sheets\Combined_CTD\combined_nuts_ctd_all_noflag_dz5.mat')

T = table();

% loop through all the combined stations
for n = 1:height(combined)
    stn = combined(n,:); % nutrient data
    % skip unnecessary stations
    if stn.combinedchar == 'Open-JS';   continue;    end  % open JS sites
    if n > 1 && any([T.StationID] == stn.StationID & [T.Year] == stn.Year); continue; end % skip stations that have been computed
   % if contains(string(stn.Location), {'to','hore'}); continue; end
    % proceed with computation
    d = contains({dts.name}, string(stn.CruiseID));
    load([ctd_root dts(d).name]); % load the relevant ctd data
    sid = strcmp([ctd.station], string(stn.ctd_name)) ;
    if sum(sid) > 1
        sid = strcmp([ctd.station], string(stn.ctd_name)) & strcmp([ctd.bottlecast], 'Y');
        if sum(sid)>1
        sid = strcmp([ctd.station], string(stn.ctd_name)) & strcmp([ctd.transect], string(stn.Transect)) ;
        end
    end
    profile = ctd(sid); % ctd data

    %% continuous variables
    % find indicies \pm 10m of the 1026 isopycnals
    [~,id_rho] = min(abs(profile.rho - 1026));
    % get outer profile
    outer = GetOutsideProf(string(stn.CruiseID));
    [~,id_rhoout] = min(abs(outer.rho - 1026));
    %{ 
    % contiuous profiles of plain anomaly
    if length(outer.depth) < length(profile.depth)
        tempa = profile.CT(1:length(outer.depth))-outer.CT';
        turba = profile.turb(1:length(outer.depth))-outer.turb';
        doa = profile.ODO_molL(1:length(outer.depth))-outer.ODO_molL';
    else
        tempa = profile.CT - outer.CT(1:length(profile.depth))';
        turba = profile.turb - outer.turb(1:length(profile.depth))';
        doa = profile.ODO_molL - outer.ODO_molL(1:length(profile.depth))';
    end
    tempa([1:id_rho-20 id_rho+20:end]) = nan; % ignore the upper 5 m 
    doa([1:id_rho-20 id_rho+20:end]) = nan; % ignore the upper 5 m 
    turba([1:id_rho-20 id_rho+20:end]) = nan; % ignore the upper 5 m 
    %}
    % compute plain anomaly around the 1026 isopycnal
    if id_rho < 20 && length(profile.depth) <  id_rho + 20  % if the profile is too shallow just take the entire watercolumn
        tempa = (mean(outer.CT(id_rhoout-20:id_rhoout+20), 'omitnan')-mean(profile.CT,'omitnan')); %(outer.CT(id_rhoout-20:id_rhoout+20), 'omitnan');
        turba = (mean(profile.turb,'omitnan')- mean(outer.turb(id_rhoout-20:id_rhoout+20), 'omitnan')); %(outer.turb(id_rhoout-20:id_rhoout+20), 'omitnan');
        doa = (mean(outer.ODO_molL(id_rhoout-20:id_rhoout+20), 'omitnan')-mean(profile.ODO_molL,'omitnan')); %(outer.ODO_molL(id_rhoout-20:id_rhoout+20), 'omitnan');
    elseif id_rho < 20 && length(profile.depth) >  id_rho + 20 
        tempa = (mean(outer.CT(id_rhoout-20:id_rhoout+20), 'omitnan')-mean(profile.CT(1:id_rho+20),'omitnan')); %(outer.CT(id_rhoout-20:id_rhoout+20), 'omitnan');
        turba = (mean(profile.turb(1:id_rho+20),'omitnan')- mean(outer.turb(id_rhoout-20:id_rhoout+20), 'omitnan')); %(outer.turb(id_rhoout-20:id_rhoout+20), 'omitnan');
        doa = (mean(outer.ODO_molL(id_rhoout-20:id_rhoout+20)-mean(profile.ODO_molL(1:id_rho+20),'omitnan'), 'omitnan')); %(outer.ODO_molL(id_rhoout-20:id_rhoout+20), 'omitnan');
    elseif length(profile) < id_rho + 20  % if 1026 depth is close to the bottom of the profile
        tempa = (mean(outer.CT(id_rhoout-20:id_rhoout+20), 'omitnan')-mean(profile.CT(id_rho-20:end),'omitnan')); %(outer.CT(id_rhoout-20:id_rhoout+20), 'omitnan');
        turba = (mean(profile.turb(id_rho-20:end),'omitnan')- mean(outer.turb(id_rhoout-20:id_rhoout+20), 'omitnan')); %(outer.turb(id_rhoout-20:id_rhoout+20), 'omitnan');
        doa = (mean(outer.ODO_molL(id_rhoout-20:id_rhoout+20), 'omitnan')-mean(profile.ODO_molL(id_rho-20:end),'omitnan')); %(outer.ODO_molL(id_rhoout-20:id_rhoout+20), 'omitnan');
    else
        tempa = (mean(outer.CT(id_rhoout-20:id_rhoout+20), 'omitnan')-mean(profile.CT(id_rho-20:id_rho+20),'omitnan')); %(outer.CT(id_rhoout-20:id_rhoout+20), 'omitnan');
        turba = (mean(profile.turb(id_rho-20:id_rho+20),'omitnan')- mean(outer.turb(id_rhoout-20:id_rhoout+20), 'omitnan')); %(outer.turb(id_rhoout-20:id_rhoout+20), 'omitnan');
        doa = (mean(outer.ODO_molL(id_rhoout-20:id_rhoout+20), 'omitnan')-mean(profile.ODO_molL(id_rho-20:id_rho+20),'omitnan')); %(outer.ODO_molL(id_rhoout-20:id_rhoout+20), 'omitnan');
    end

    %% nitrate anomaly
    % all depths at stn
    id_stn = combined.StationID == stn.StationID & combined.Year == stn.Year;
    comb_stn = sortrows(combined(id_stn,:), 'Depth_m'); 
    % get a "mean" nutrient profile from stations in open jones sound, use
    % broad depth bins
   
    %id_surf = find(comb_stn.Depth_m <= 30); id_deep = find(comb_stn.Depth_m > 30); dno3a = -(mean(comb_stn.NO3_uM(id_deep),'omitnan') - mean(comb_stn.NO3_uM(id_surf),'omitnan'));
    dno3 = diff(comb_stn.NO3_uM); id_surf = find(comb_stn.Depth_m <= 30); dno3a = mean(dno3(id_surf(1:end-1)), 'omitnan');
    no3 = mean(comb_stn.NO3_uM(id_surf)); nh3 = mean(comb_stn.NH3_uM(id_surf));

    %% save into table
    T1 = table();

    T1.tempa = tempa; T1.turba = turba; T1.doa = doa; T1.dz = stn.dz;
   % T1.depth = profile.depth;
    T1.StationID = stn.StationID; T1.Year = stn.Year; T1.Transect = stn.Transect; T1.Location = stn.Location;
    T1.dno3 = dno3a; T1.fwsource = stn.FWSource; T1.chla = stn.chla_c;
    T1.no3 = no3; T1.turb = stn.turb_c; T1.dist = stn.dist; T1.NH3 = nh3;
    %T1.mn_surf = mean(comb_stn.NO3_uM(id_surf),'omitnan'); T1.mn_deep = mean(comb_stn.NO3_uM(id_deep),'omitnan');
    
    if n == 1
        T = T1;
    elseif n >1
        T = [T; T1];
    end
end

% scale to [-1 1] scale
vars = {'doa';'tempa';'turba';'dz';'dno3'};
for vv = 1:length(vars)
   T.([vars{vv} 'N']) = normalize(T.(vars{vv}), 'medianiqr');
   %T.([vars{vv} 'N']) = normalize(T.(vars{vv}), 'zscore', 'robust');
   %T.([vars{vv} 'N']) = -1 + 2*(T.(vars{vv}) - min(T.(vars{vv})))./(max(T.(vars{vv}))- min(T.(vars{vv})));
end

T.idx_N =  T.doaN + +T.tempaN + T.dzN + T.dno3N + T.turbaN;
T.idx =  T.doaN + +T.tempaN + T.dzN + T.turbaN;

%
save('C:\Users\claire\Dropbox (Bhatia Lab)\Bhatia Lab Team Folder\Projects\Devon\Multi-year_analyses\Analysis\nutrients\master_sheets\combined_CTD\UpwellingIndex.mat', 'T');
%writetable(T, ['C:\Users\claire\Dropbox (Bhatia Lab)\Bhatia Lab Team Folder\Projects\Devon\Multi-year_analyses\Analysis\nutrients\master_sheets\combined_CTD\UpwellingIndex.csv'])
figure;
bar(1:height(T), T.idx); hold on;

T2 = T(~isnan(T.idx),:);
figure;
boxchart(T2.Location, T2.idx, 'GroupByColor',T2.Year)