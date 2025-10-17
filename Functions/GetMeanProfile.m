function [prof] = GetMeanProfile(data, var)
%%
% Calculate mean profile for given structure
%  
% data Mx1 profile
% 
%%%
% get the max length to set the empty profile length
max_len = 0;
if isempty(data) 
    prof.(var) = []; 
    prof.([var, 'std']) = []; 
    prof.depth = []; 
    return
end 
for n = 1:length(data)
    len = length(data(n).depth);
    if len > max_len
        max_len = len; 
        dep = data(n).depth;
    elseif len == 0
        prof = [];
        return 
    end
end
clear n len

% Set up empty matrix to populate
profarr = ones(length(data), max_len)*nan;
for n = 1:length(data)
    if strcmp(var, 'fwt')
        data_fw = CalcFWT(data(n), 34.8, 100, 'profile');
        id = length(data_fw);
        profarr(n,1:id) = data_fw;
    elseif strcmp(var, 'ce')
        data_ce = CalcCEProf(data(n), 100);
        id = length(data_ce);
        profarr(n,1:id) = data_ce;
    else
        id = length(data(n).(var));
        profarr(n,1:id) = data(n).(var);
    end
end

% save in structure
prof.(var) = mean(profarr,1,'omitnan');
prof.([var, 'std']) = std(profarr,1,'omitnan');
prof.depth = dep;
end