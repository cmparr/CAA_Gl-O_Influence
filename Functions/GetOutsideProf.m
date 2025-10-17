function [cast ctdf] = GetOutsideProf(cruiseID, varargin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Get a profile from outside coastline to use as background to look for
%   coastline anomalies. Gives a mean profile of all 'outer' profiles from
%   samples taken in 2019-2022.
%   
%   INPUT:
%       data    MxN Struct
%               data set containing profile
%
%   OUTPUT:
%       cast    1XN Struct
%               selected background profile
%       ctdf    MXN Struct
%               structure of all data going into the mean cast with M
%               profiles. 
%
%   A space saving maneuver 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ids = {'V1907'; 'V2007'; 'V2108'; 'A2108'; 'V2208'; 'A2208'};
idx = strcmp(ids, cruiseID);
yr_ids = {'2019'; '2020'; '2021'; '2021'; '2022'; '2022'};
yr = yr_ids{idx}; id_yr = find(strcmp(yr, yr_ids));
% No bottom depth for small boat campaigns, use given profiles
A2108_stnuse = {'JS1'; 'JS2';'JS8'; 'JS52'};
A2208_stnuse = {'JSOG1B'};
variables = {'press';'temp';'cond';'sal';'turb'; 'chla'; 'ODO_molL'; 'ODO_perc'; 'depth'; 'PAR';'SA';'CT';'pt';'sigma0';'rho';'N2';'N2press';'CE'};
%stns = {'VIO33';'JSD';'JS12';'JS12';'JS4'; 'JS4'}; 
%
%cruiseID = ids_use{idx};
%% Get file, stn name based on cruiseID
dts = dir('C:\Users\claire\Dropbox (Bhatia Lab)\Bhatia Lab Team Folder\Projects\Devon\Multi-year_analyses\Annual_CTD_Processed\Annual_CTD_matfiles\Claire_PROCESSED\*.mat');
filename = dts(contains({dts.name}, cruiseID)).name;
load(['C:\Users\claire\Dropbox (Bhatia Lab)\Bhatia Lab Team Folder\Projects\Devon\Multi-year_analyses\Annual_CTD_Processed\Annual_CTD_matfiles\Claire_PROCESSED\' filename]);

ctdf = ctd([ctd.bottom_depth] >= 950.); % get an empty structure like each ctd 
%% get the profiles for all casts deeper than 250m depth
for ii = 1:length(id_yr)
    cruiseid = ids(id_yr(ii));
    filename = dts(contains({dts.name}, cruiseid)).name;
    load(['C:\Users\claire\Dropbox (Bhatia Lab)\Bhatia Lab Team Folder\Projects\Devon\Multi-year_analyses\Annual_CTD_Processed\Annual_CTD_matfiles\Claire_PROCESSED\' filename]);
    if strcmp(cruiseid, {'A2108'})
        for ss = 1:length(A2108_stnuse)
            ctd1 = ctd(strcmp([ctd.station], A2108_stnuse{ss}));
            ctdf = [ctd1 ctdf];
        end
    elseif strcmp(cruiseid, {'A2208'})
       ctd1 = ctd(strcmp([ctd.station], A2208_stnuse));
       ctdf = [ctd1 ctdf];
    elseif  strcmp(cruiseid, {'V2208'})
        ctd1 = ctd([ctd.bottom_depth] >= 250. & ~contains([ctd.transect], 'Talbot') & ~contains([ctd.transect], 'Grise')& ~strcmp([ctd.station], '2021INC') & ~contains([ctd.transect],'Belcher') );
        ctd1 = ctd1(~strcmp([ctd1.station], 'OG') & ~strcmp([ctd1.station], 'OG_A') & ~contains([ctd1.bottlecast],'Y') );
        ctdf = [ctd1 ctdf];
    elseif strcmp(cruiseid, {'V2208'})&strcmp(varargin, 'Talbot')
        ctd1 = ctd([ctd.bottom_depth] >= 250. & strcmp([ctd.station], 'TA1'));
        ctdf = [ctd1 ctdf];
    else
        ctd1 = ctd([ctd.bottom_depth] >= 250. & ~contains([ctd.transect], 'Talbot') & ~contains([ctd.transect], 'Grise')& ~strcmp([ctd.station], '2021INC') & ~contains([ctd.transect],'Belcher') );
        ctdf = [ctd1 ctdf];
    end
end

%% take average for watercolumn variables
cast = struct();
% metadata
cast.cruiseids = ids{id_yr};
cast.year = yr;

for vv = 1:length(variables)
    meanstruct = GetMeanProfile(ctdf, variables{vv});
    cast.(variables{vv}) = meanstruct.(variables{vv});
    cast.([variables{vv}, 'std']) = meanstruct.([variables{vv}, 'std']);
    if strcmp(variables{vv}, 'depth')
        cast.(variables{vv}) = meanstruct.(variables{vv})';
    end
end
end
