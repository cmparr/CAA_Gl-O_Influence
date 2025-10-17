function [dz] = Getdz(station, location_name, cruiseID)
% Calculate the upwelling index for the given station and transect. Is
% physically the difference in depth of the 1026 isopycnal at the outermost
% profile and the depth of the 1026 isopycnal at the station. 
%
%   INPUT:      station
%               station name from ctd (note not from nutrients data - they are different)               
%
%   OUTPUT:     dz = station depth of upwelling


%% get the full ctd profile based on cruiseid
dts = dir('C:\Users\claire\Dropbox (Bhatia Lab)\Bhatia Lab Team Folder\Projects\Devon\Multi-year_analyses\Annual_CTD_Processed\Annual_CTD_matfiles\Andrew_PROCESSED\*.mat');
filename = dts(contains({dts.name}, cruiseID)).name;
load(['C:\Users\claire\Dropbox (Bhatia Lab)\Bhatia Lab Team Folder\Projects\Devon\Multi-year_analyses\Annual_CTD_Processed\Annual_CTD_matfiles\Andrew_PROCESSED\' filename]);

%% Isolate the transect in the CTD file, sort based on distance from coastline
transect = struct();

idx = contains(string([ctd.location]), location_name).*contains(string([ctd.transect]),'out'); 
if contains(location_name,'Jakeman practice')
    idx = contains(string([ctd.transect]),'Jakeman Glacier - practice'); 
elseif contains(location_name, 'Sydkap transect') & cruiseID == 'V1907'
    idx = contains(string([ctd.location]), 'South Cape').*contains(string([ctd.transect]),'out'); 
elseif contains(location_name, 'Sydkap Sill transect') & cruiseID == 'A2108'
    idx = contains(string([ctd.location]), 'Sydkap').*contains(string([ctd.transect]),'out'); 
elseif station == 'VIO16' | station == 'VIO14'
    idx = contains(string([ctd.transect]),'Belcher Shore - out'); 
elseif strcmp(location_name,'Sverdrup Glacier-repeat')
    idx = contains(string([ctd.location]), 'Sverdrup Glacier').*contains(string([ctd.transect]),'out'); 
end
if sum(idx) == 0
    dz = nan; id_stn = nan; id_iso = nan; dz_transect = nan;
end

if length(transect) > 1
   transect = [transect ctd(logical(idx))];
else
   transect = ctd(logical(idx));
end
if isempty(transect);   return;  end

[dists, dist_idcs] = find_distance(transect); transect = transect(dist_idcs);
%{
if min(dists) < 5.0  % previous filtering for stations within 5km from
coastline
    in_transect = transect(dists <= 5.0);
else
    in_transect = transect(1); % if no near-glacier measurements take inner-most station
end
%}

%% Get the outermost profile & isopycnal depth
%{
if contains(transect_name, 'Sverdrup')&& contains(transect_name, 'out')
    id_out = length(in_transect) - 1; % use site almost closest to inside of moraine at Sverdrup
elseif contains(transect_name, 'Jakeman')&& contains(transect_name, 'out')&& ~strcmp(cruiseID, 'V1907')
    [~, id_out] = min(abs(dists - 3.0)); % use site almost closest to inside of moraine at Jakeman
else
    id_out = length(in_transect); 
end
%}

id_out = length(transect); 

% index of 1026 isopycal for outer profile 
[~, id_26o] = min(abs(transect(id_out).rho - 1026));

%% Find station profile
id_stn = strcmp([transect.station], station); 
if sum(id_stn)==0; 
    id_stn = strcmp([ctd.station], station);  
    if sum(id_stn) > 1
        id_stn = strcmp([ctd.station], station)&strcmp([ctd.bottlecast], 'Y');
    end
    stn = ctd(id_stn);
elseif sum(id_stn)==1; stn = transect(id_stn); 
elseif sum(id_stn)>1 & all([transect.cruiseID] == 'V2108');
    id_stn = strcmp([transect.station], station)& strcmp([transect.transect], 'Talbot Inlet - out');
    stn = transect(id_stn); 
elseif sum(id_stn)>1 & all([transect.cruiseID] == 'V2208');
    id_stn = strcmp([transect.station], station)& strcmp([transect.bottlecast], 'Y');
    stn = transect(id_stn); 
end


[~, id_26s] = min(abs(stn.rho - 1026));

dz = transect(id_out).depth(id_26o) - stn.depth(id_26s);

end