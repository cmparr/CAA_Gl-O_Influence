%%% Script to run PCA analysis in Thesis/Paper. Includes data cleaning,
%%% analysis and plotting. Figures 4 and 7.
%%% Written by Claire Parrott (cparrott@eoas.ubc.ca), last modified 13/01/2025

close all
clear

direct = dir('../Data/*CTD_Data_50cm_binned.mat');

%% set some figure properties
set(0, 'defaultAxesFontSize', 15);

%% isolate just the numeric data from the table; remove bad turbidity measurements
conditions = {' Full Column'; ' Upper 30m'; ' Below 30m'};
fig_labels = {{'';'$(a)$';'$(c)$'};{'';'$(b)$';'$(d)$'}};

subset_data =  {'allsites'; 'glaciers'}; % which subset of data to use

for c = 2 %2:3 % upper or lower wcol
cond = conditions{c}; % filtering condition applied to input data
%
for d = length(direct)
    cd('./Data');
    load(direct(d).name);

ctypes = {'space';'time';'continuum'};
for y = 1:2  % spatial or year or continuum
    ctype = ctypes{y};
    %fig_num = fig_labels{y};

if d <5 
cond = strcat(' ', direct(d).name(19:22), cond); % filtering condition applied to input data
end 

combined = combined((combined.dist < 20.) | (combined.Location == 'Jones Sound'), :);
% omit sites
idx = find(~strcmp(string([combined.Transect]), 'Sydkap Glacier - across') & ~contains(string([combined.Location]), {'Talbot'; 'Harbour';'Terry';'to';'Shore'}));
combined = combined(idx,:);

for s = 1 % 1 for all sites, 2 for glaciers
    subset = subset_data{s};
%
if s == 2
idx = find(~contains(string([combined.Location]), {'Grise'; 'Fram'; 'True';'Terry'}));
combined = combined(idx,:);
y = 3;
end
%% select variables
% set variable names
variables =  {'CT', 'sa', 'CE30', 'DO_Lmol','chla', 'turb', 'NO3_uM','PO4_uM', 'SiO4_uM'};
variable_names = {'Temperature', 'Salinity','CR', 'Dissolved O$_2$','Chlorophyll-$\alpha$', 'Turbidity', 'NO$_3$', 'PO$_4$', 'SiO$_4$'};
if s == 2
variables = {'CT', 'DO_Lmol','chla', 'turb', 'NO3_uM', 'PO4_uM', 'SiO4_uM', 'dz'};
variable_names = {'Temperature','Diss. O$_2$', 'Chlorophyll-$\alpha$', 'Turbidity', 'NO$_3$', 'PO$_4$', 'SiO$_4$', '$\Delta$z'};
end

% depth partitioning, create matrix to be used in PCA
% omit nan values
idd = ~any([isnan(combined.NO3_uM), isnan(combined.SiO4_uM), isnan(combined.PO4_uM)]');
combined = combined(idd,:);
idd = any([~isnan(combined.NO3_uM), ~isnan(combined.SiO4_uM), ~isnan(combined.PO4_uM)]');
combined = combined(idd,:);

if contains(cond, ' Upper 30m')
    combined = combined(combined.Depth_m <= 30., :); % filter out anything at depth > 30m 
    data = [];
    for vv = 1:length(variables)
        data = [data combined.(variables{vv})];
    end
    disp('upper 30')
elseif contains(cond, ' Below 30m')    
    combined = combined(combined.Depth_m > 30. & combined.Depth_m <= 100., :); % filter out anything at depth < 30m 
    data = [];
    for vv = 1:length(variables)
        data = [data combined.(variables{vv})];
    end
    disp('below 30')
else
    data = [];
    for vv = 1:length(variables)
        data = [data combined.(variables{vv})];
    end
end

%% perform PCA
rowidx= any(isnan(data'));
data = data(~rowidx,:);
combined = combined(~rowidx',:);
[ei, pcs, eigv, var, ~] = CenteredPCA(data, variables, 'plot');
fig_var = figure(1); fig_eig = figure(2);
close 3; % dont need the scatter plot making a better one
p1 = 1; p2 = 2; % integers for each pc axis used

%{
%% Test significance of groupings, takes a while better not to run if you know the results
[pg, site_groups] = perMANOVA(data, combined.combinedchar)
if s == 2
[pgl, site_glaciers] = perMANOVA(data, combined.Location)
end
[py, yr_groups] = perMANOVA(data, combined.Year)
%}
%% plot PC values in PC1 PC2 space 
fig4 = figure(4); % axes for plotting custom scatter plot
fig4.Position= [50 50 700 700];
hold on;

%% Format eigenvectors/background
% plot eigenvectors
% Set quiver of eigenvectors as base
EI = 8.5*ei; 
off = 0.1.*[EI(:,p1)./abs(EI(:,p1)), EI(:,p2)./abs(EI(:,p2))]';
%
% offset eigenvector labels
if s == 1
%  vars:  {'CT', 'SA', 'CR','DO','Chla', 'Turb', 'NO3','PO4', 'SiO4'};  ;
    if c == 2
     off2 = [-0.9 -1 -0.1 -0.5 -1 -0.6 -0.3 -0.3 -0.3; % upper 30 BGC, values selected by hand
            0.0 0.0 -0.3 -0.2 -0.2 -0.1 0.1 0.1 -0.1];
    elseif c == 3
    off2 = [-0.4 -0.4 -0.4 -0.4 -0.5  -0.7  -0.2 -0.2 -0.2; %below 30m BGC
            0.0  0.2  -0.5  0.1  0.1  -0.1  0.0  0.0  0.0];
    else
    end
%  vars:  {'CT', 'DO_Lmol','chla', 'turb', 'NO3_uM', 'PO4_uM', 'SiO4_uM', 'dz'};
elseif s ==2
    if c == 2
    off2 = [-1.1 0.0 -1.0 -0.4 -0.4 -0.4 -0.4 -0.3 ; %upper 30m glaciers
            0.0 -0.3 -0.4 -0.3  0.2  0.1 -0.2 -0.4];
    elseif c == 3
    off2 = [-0.4 -0.4 -0.2 -0.7 -0.0 0.0 ; %below 30m glaciers
            -0.0 -0.3  0.0  0.1 -0.2 -0.2 ];
    else 
    end
end
off = off + off2;
%}
% guide lines, plot eigenvector labels
xline(0, 'HandleVisibility','off'); yline(0, 'HandleVisibility','off'); hold on;
q = quiver(zeros(size(EI(:,p1))), zeros(size(EI(:,p1))), EI(:,p1), EI(:,p2), 'k', 'Linewidth', 1, 'HandleVisibility','off'); 
text(EI(:,p1)+off(1,:)', EI(:,p2)+off(2,:)', variable_names, 'FontSize', 13, 'Interpreter','latex');

%% set color matrix to plot by different features
colors = zeros([3 size(combined,1)]); shapes = zeros([size(combined,1) 1]);
% color by type of site
disp(direct(d).name(19:22))
disp(size(combined))
twstrong = contains(string([combined.Location]), {'Belcher'; 'Sverdrup'; 'Jakeman'}); % glacierized & open
twweak   = contains(string([combined.Location]), {'Sydkap','South Cape'}); % glacierized & fjorded
riverine = contains(string([combined.Location]), {'True'}); % riverine, open
dry      = contains(string([combined.Location]), {'Grise'; 'Fram'}); % riverine, fjord
open_sites=contains(string([combined.Location]), {'Jones Sound';'OG'}); 

d=5; 
%% set colors 
% shape by year of observations
y2019 = combined.Year == 2019; 
y2020 = combined.Year == 2020; 
y2021 = combined.Year == 2021; 
y2022 = combined.Year == 2022; 
markers = {'s'; 'o'; '^'; 'pentagram'};
d=5; 
%} 

if y == 1 % spatial colors
    colors(:,twstrong) = repmat([0; 0.447; 0.8], [1 sum(twstrong)]); 
    colors(:, twweak)  = repmat([0.850; 0.325; 0.098], [1 sum(twweak)]);   
    colors(:, riverine) = repmat([0.466; 0.647; 0.188], [1 sum(riverine)]); 
    colors(:, dry)     = repmat([0.929; 0.694; 0.125], [1 sum(dry)]);  
    colors(:, open_sites)     = repmat([1; 1; 1]*0.6, [1 sum(open_sites)]);  
    %

    %% legend for colors
    id_cols = [find(twstrong,1) find(twweak,1) find(riverine,1) find(dry,1) find(open_sites,1)]; id_cols = id_cols(id_cols > 0);
    
    for ii = 1:length(id_cols)
    scatter(pcs(id_cols(ii),p1), pcs(id_cols(ii),p2), 1, colors(:, id_cols(ii))', 'filled', 'marker', 'o', 'HandleVisibility','on'); hold on;
    end
    %{    % multi-year
    for yy = 1:4
    scatter(-20, -20, 4, [1 1 1]*0.7, 'filled', 'marker', markers{yy}, 'HandleVisibility','on','MarkerFaceAlpha',1, 'MarkerEdgeColor', 'k'); hold on;
    end
    %
    %leg_fill = {'Glacierized-Open'; 'Glacierized-Fjorded'; 'Riverine-Open'; 'Riverine-Fjorded'};
    leg_fill = {'Glacierized-Open'; 'Glacierized-Fjorded'; 'Riverine-Open'; 'Riverine-Fjorded'; 'Open Jones Sound'; '2019';'2020';'2021';'2022'};
elseif y == 2  %% years colors
    colors = cmocean('haline', 4)';
    colors(:, 4) = brighten(colors(:,4), -0.4);
    id_cols = [find(y2019,1) find(y2020,1) find(2021,1) find(2021,1)]; id_cols = id_cols(id_cols > 0);
    for ii = 1:length(id_cols)
       scatter(pcs(id_cols(ii),p1), pcs(id_cols(ii),p2), 1, colors(:, ii)', 'filled', 'marker', 'o', 'HandleVisibility','on'); hold on;
    end
    %
    leg_fill = num2str([2019:2022]');
elseif y == 3 % continuum
    %% glacier continuum
    colors = zeros([3 size(combined,1)]); shapes = zeros([size(combined,1) 1]);
    twstrong =  contains(string([combined.Location]), {'Belcher'}); % 
    twmid1    = contains(string([combined.Location]), {'Sverdrup'});
    twmid2    = contains(string([combined.Location]), {'Sydkap';'South Cape'});
    twweak   =  contains(string([combined.Location]), {'Jakeman'}); % 
    riverine =  contains(string([combined.Location]), {'True'; 'Grise'; 'Fram'}); % riverine, open
    open_sites =contains(string([combined.Location]), {'Jones Sound'});
    %
    cmap = cmocean('matter', 4);
    % 
    colors(:,twstrong) = repmat(cmap(1,:)', [1 sum(twstrong)]); 
    colors(:, twmid1)  = repmat(cmap(3,:)', [1 sum(twmid1)]);   
    colors(:, twmid2)  = repmat(cmap(2,:)', [1 sum(twmid2)]);   
    colors(:, twweak)  = repmat(cmap(4,:)',[1 sum(twweak)]);   
    colors(:, riverine)     = repmat([0.466; 0.647; 0.188], [1 sum(riverine)]);  
    colors(:, open_sites)     = repmat([1; 1; 1]*0.6, [1 sum(open_sites)]);

   %% legend for colors
    id_cols = [find(twstrong,1) find(twmid2,1) find(twmid1,1) find(twweak,1) find(riverine,1) find(open_sites,1)]; id_cols = id_cols(id_cols > 0);
    
    %% plot for legend values
    for ii = 1:length(id_cols)
    scatter(pcs(id_cols(ii),p1)+20, pcs(id_cols(ii),p2)+20, 1, colors(:, id_cols(ii))', 'filled', 'marker', 'o', 'HandleVisibility','on', 'MarkerEdgeColor', 'k'); hold on;
    end
    for yy = 1:4
    scatter(-20, -20, 4, [1 1 1]*0.7, 'filled', 'marker', markers{yy}, 'HandleVisibility','on','MarkerFaceAlpha',1, 'MarkerEdgeColor', 'k'); hold on;
    end
    %
    leg_fill = {'Belcher'; 'Sydkap'; 'Sverdrup'; 'Jakeman'; 'Open Jones Sound';'2019';'2020';'2021';'2022'};
end
% set legend
legend(leg_fill, 'Location','northeast','NumColumns',2);

%% plot data
if d == 5
    yrs = [y2019 y2020 y2021 y2022]; 
    for yy = 1:4
        data = [pcs(yrs(:,yy),p1)  pcs(yrs(:,yy),p2)];
        if y == 1 || y == 3
        scatter(pcs(yrs(:,yy),p1), pcs(yrs(:,yy),p2), 120, colors(:,yrs(:,yy))', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7, 'HandleVisibility', 'off', 'marker', markers{yy}); hold on; 
        elseif y == 2
        scatter(pcs(yrs(:,yy),p1), pcs(yrs(:,yy),p2), 120, colors(:,yy)', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7, 'HandleVisibility', 'off', 'marker',  'o'); hold on; %plot yrs in scatter
        end
        %[~, avg, lens] = PlotErrorEllipse(data,'color', colors(:,yy)');  % color ellipses with year
        %avgs(yy,:) = avg; 
    end
else
scatter(pcs(:,p1), pcs(:,p2), combined.dist*7, colors', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.8, 'marker', 'o', 'HandleVisibility','off'); hold on;
end
%{
cat_idx = [twstrong twweak riverine dry];
cat_name = {'twstrong'; 'twweak'; 'riverine'; 'dry'};
for cc = 1:size(cat_idx,2) % color ellipses by site classif
cidx = cat_idx(:,cc);
    [~, avg, lens] = PlotErrorEllipse(pcs(cidx,p1:p2), 'color', colors(:,id_cols(cc)), 'save','cond', [cond cat_name{cc}]); 
end
%}

%% figure formatting
daspect([1 1 1]);
xlabel(strcat('PC', num2str(p1),  '  (', num2str(var(p1)*100, '%.1f'), '%)'));
ylabel(strcat('PC', num2str(p2),  '  (', num2str(var(p2)*100, '%.1f'), '%)'));
xlim([-4.5 7]); ylim([-3.5 6.5]);
grid();

%% saving
 rootfol = '../figures/';
%
conds = string(regexprep(cond, ' ', ''));
ax = gca;
saveas(fig4,strcat(rootfol,'PC', num2str(p1), 'vsPC', num2str(p2), '-scatter-', ctype,'-', subset, '-', conds,'-opensites.png'));


close all
cond =  conditions{c};
%}
end
end
end
end



