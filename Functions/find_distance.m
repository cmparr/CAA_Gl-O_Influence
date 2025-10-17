function [distances, sorted_idcs] = find_distance(data, varargin)
% FIND_DISTANCE Distance from one coordinate to common point
% ========================================================================
% USAGE: [distances] = find_distance(data)
%
% DESCRIPTION: Find the distance between a set of coordinates and
% a common reference point for each transect. Uses sw_dist, finds Plane
% Sailing method between two points
%
% INPUT: data (data NxM) 
%        varargin: 'full' - For the full across Jones Sound transect, set
%        1st point to Grise Fiord
%              'talbot-full' - use the to-glacier distance for transect for
%              Talbot transect.
%
% OUTPUT: distances, sorted_idcs    (double 1XM, indicies 1XM) 
%           Distances of each point from common location with sorted index
%%%%%%

load('C:\Users\claire\Dropbox (Bhatia Lab)\Bhatia Lab Team Folder\Projects\Devon\Multi-year_analyses\Annual_CTD_Processed\Annual_CTD_mfiles\Claire_mfiles\Functions\commonglacierloc.mat')
cgl = CommonGlacierLocation;
if istable(data)
    len = height(data);
    trans = 'Transect';
    la = 'Lat_dd_ddddN_'; lo = 'Lon_dd_ddddE_';
else
    len = length(data);
    trans = 'transect';
    la = 'lat'; lo = 'lon';
end
distances = zeros(len,1);

for ii=1:len
    if istable(data)
        data1 = data(ii,:);
    else
        data1 = data(ii);
    end
    tname = string(data1.(trans));
    if strcmp(tname,'Central Jones Sound')|strcmp(tname,'Western Entrance Jones Sound')|strcmp(tname,'Eastern Jones Sound Deep Water')
        tname = 'Jones Sound';
    end
    if any(strcmp(varargin, 'full'))
        tname = 'Grise - full';
    elseif any(strcmp(varargin, 'talbot-full'))
        tname = 'Talbot-True';
    end
  
    if contains(tname, cgl.Location)
            if strcmp(tname, 'Grise - full')
                cidx = find(contains(cgl.Location, extractBetween(tname, 1, 4)));
                cidx = cidx(2);
            elseif contains(tname, 'Talbot Inlet')
                cidx = find(contains(cgl.Location, extractBetween(tname, 1, 4)));
                cidx = cidx(2);
            else
                cidx = find(contains(cgl.Location, extractBetween(tname, 1, 4)),1);
            end
        reflat = cgl.Lat(cidx);
        reflon = cgl.Lon(cidx);
    else
        reflat = data(1).(la);
        reflon = data(1).(lo);
    end
        distances(ii) =sw_dist([reflat data1.(la)],[reflon data1.(lo)],'km'); % dist of each station from first station in set
end   

distances = round(distances,3);
[~, sorted_idcs] = sort(distances);
end








