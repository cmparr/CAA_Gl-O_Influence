
%Import Standardized Metadata CSV files
%Important to convert xls files to csv or there will be issues with datetime
%formats - this might not be an issue anymore with new code (should revisit
%- maybe can import direct from xls?)

%Run this m-file first before running 'A2_merge_metadata_rsk_raw.m'
%Then run 
%Created by Andrew Hamilton (akhamilt@ualberta.ca)
%Created Dec 7, 2022
%Last modified Dec 7, 2022

%Read metadata sheet and find breakline that delineates header from fieldnames and data rows
clearvars
rootfol = "/Users/claire/Dropbox (Bhatia Lab)/Bhatia Lab Team Folder/Projects/Devon/Multi-year_analyses/Annual_Logs_Standardized/Annual_CTD_Logs_Standardized_csv/";
cd(rootfol);
d = dir('*.csv');
for n = 1:length(d)
    csvfile = d(n).name
    headerbreak = '###*'; %Breakline contains series of # characters
    fid=fopen(csvfile,'r');
    for nn = 1:100; %Read up to 100 lines in csv file to find break characters
        strline=fgetl(fid);
        if regexp(strline,headerbreak);
            breakline = nn;
            break %Stops loop when ### characters reached
        end
    end
fclose(fid);
clearvars -except breakline d n csvfile rootfol

%Read in metadata table from csv
opts = detectImportOptions(csvfile);
opts.Delimiter = ",";
opts.VariableNamingRule = 'preserve';
opts.VariableNamesLine = breakline+1;
opts.DataLines = breakline+2;
metadata = readtable(csvfile, opts);

%Remove special charcters from variable names
varnames  = metadata.Properties.VariableNames;
%varnames = regexprep(string(varnames), {'_'}, {''});
varnames = regexprep(string(varnames), {'[().:_\- /]+'}, {''});
metadata.Properties.VariableNames = varnames;

%Change metadata variable types to strings, double only
%ID which fields we want as doubles (all the rest will be converted to strings)
%This may not work if fieldnames changes - but won't loose data, just might
%not be in preferred data type (easily converted as needed)
doublevars = {'CastNo', 'LatddddddN',	'LonddddddE',	'BottomDepthm',	'MaxCTDDepthm',...
    'Bottle1Depthm' , 'Bottle2Depthm'	, 'Bottle3Depthm', 'Bottle4Depthm' , 'Bottle5Depth_m' , ...
    'Bottle6Depthm' , 'Bottle7Depthm' , 'Bottle8Depthm' ,'Bottle9Depthm' , 'Bottle10Depth_m' ,...
    'RSKFullProfileno'	,	'Icethicknesscm'	, 'Icefreeboardcm',	'Snowthicknesscm'};
logdouble = ismember (varnames, doublevars);
inddouble = find(logdouble);
logstring = ~logdouble; %Any field that is not supposed to be a double will be a string
indstring = find(logstring);
%Convert number fields from char to doubles
for n = 1:length(inddouble);
        if iscell(metadata.(varnames{inddouble(n)}));
        metadata.(varnames{inddouble(n)}) = str2double(metadata.(varnames{inddouble(n)}));
    end
end
%Convert char and datetimes and durations to strings
for n = 1:length(indstring);
        if iscell(metadata.(varnames{indstring(n)})) | isa(metadata.(varnames{indstring(n)}), 'double') | isa(metadata.(varnames{indstring(n)}), 'datetime') | isa(metadata.(varnames{indstring(n)}), 'duration');
        metadata.(varnames{indstring(n)}) = string(metadata.(varnames{indstring(n)}));
    end
end

filename = [csvfile(1:5) '_CTD_Metadata.mat']
cd('/Users/claire/Dropbox (Bhatia Lab)/Bhatia Lab Team Folder/Projects/Devon/Multi-year_analyses/Annual_Logs_Standardized/Annual_CTD_Logs_Standardized_mat/');
save(filename, 'metadata')
cd(rootfol)
clearvars -except d n rootfol
end
