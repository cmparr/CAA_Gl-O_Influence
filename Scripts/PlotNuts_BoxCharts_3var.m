%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot box and whisker plot for macronutrient concentrations including
% t-test results. Figure 5.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
close all

%% set some figure properties
set(0, 'defaultAxesFontSize', 13);

% load data
load('../Data/combined_nuts_ctd_all_noflag_dz5.mat');

LRce = struct(); LRfwt = struct(); ll = 1;
years = 2019:1:2022;
% Select desired sites
%combined = combined((combined.dist < 20. & contains(combined.Location, {'Belcher'; 'Sverdrup'; 'South Cape';'Jakeman'}) & ~contains(combined.Location, {'Shore';'Talbot';'Point';'Terry';'OG'})), :);
combined = combined((combined.dist < 20. & ~contains(combined.Location, {'Jones Sound'; 'Harbour'; 'Shore';'Talbot';'Point';'Terry';'OG'})), :);

% create categorical varibales to help plotting
combined.FWSourceC = categorical(combined.FWSource); combined.CoastC = categorical(combined.Coastline); combined.CombcharC = categorical(combined.combinedcharshort);
combined.depth_bin = repmat(">100m", [height(combined) 1]); combined(combined.Depth_m <= 30.,:).depth_bin = repmat("Upper 30m", [height(combined(combined.Depth_m <= 30.,:)) 1]);  combined(combined.Depth_m > 30.&combined.Depth_m <= 100.,:).depth_bin = repmat("30-100m", [height(combined(combined.Depth_m > 30.& combined.Depth_m <= 100.,:)) 1]); 
combined.depth_bin = categorical(combined.depth_bin);
combined.CEd = combined.CE; combined(combined.Depth_m <= 30.,:).CEd = combined(combined.Depth_m <= 30.,:).CE30;  combined(combined.Depth_m > 30.&combined.Depth_m <= 100.,:).CEd = combined(combined.Depth_m > 30.& combined.Depth_m <= 100.,:).CE100; 
ch_name = char(combined.Location); 
combined.name_short = categorical(cellstr(ch_name(:,1:4)));

% figure elemetns
variables = {'NO3_uM','PO4_uM', 'SiO4_uM'};
labels = {'NO_3 Concentration (\mu M)'; 'PO_4 Concentration (\mu M)'; 'SiO_4 Concentration (\mu M)'};
titles = {'Nitrate'; 'Phosphate'; 'Silicate'};
fig_nums = {'a';'b';'c';'d';'e';'f'};
lims = [[0 12]; [0 1.2]; [0 20]];
ri_means = [0.184, 0.267, 1.573];

%% ok now make some box charts
ax = figure(1); ax.Position = [0 0 1400 600];
tiledlayout(2,3,'TileSpacing','compact', 'Padding','compact');
%combined = combined(combined.Depth_m > 30.&combined.Depth_m <= 100.,:);
xx = 0;
stars = {'***';'****';'****';'';'';''};
for yy = 1:6
    t = nexttile(yy);
    % loops based on tile so correct variale is plotted
    if yy <= 3
        comb = combined(combined.Depth_m < 30.,:);
        xx = xx +1;
    elseif yy == 4
        comb = combined(combined.Depth_m > 30. & combined.Depth_m < 100.,:);
        xx = 1;
    elseif yy > 4
        comb = combined(combined.Depth_m > 30. & combined.Depth_m < 100.,:);
        xx = xx +1;
    end

    boxchart(comb.FWSource, comb.(variables{xx}), 'BoxEdgeColor', 'k', 'BoxFaceAlpha',0.5, 'MarkerStyle', '.', 'MarkerSize', 20); hold on;
    scatter([1 2], [mean(comb(comb.FWSource=='Glacier',:).(variables{xx}), 'omitnan') mean(comb(comb.FWSource=='River',:).(variables{xx}))], 50, 'k', 'filled', 'o', 'HandleVisibility','off');
    % t-test
    [h,p,ci,stats] = ttest2_cp(comb(comb.FWSource=='Glacier',:).(variables{xx}),comb(comb.FWSource=='River',:).(variables{xx}));
    xline(1.5, ':','LineWidth', 1, 'HandleVisibility','off');
    ylim(lims(xx,:)); 
    ylabel(labels{xx}); 
    xaxis=get(gca,'XAxis');
    xaxis.TickLabels = {'Glacierized'; 'Riverine'};
    if p < 0.02; x = 2.; rr = '%.2d'; else x = 2.2; rr = '%0.2f'; end
    text(x, lims(xx,2)-lims(xx,2)/20,['p-value: ' num2str(p,rr)], fontsize = 12);
    text(x, lims(xx,2)-lims(xx,2)/30,stars{yy}, FontSize=20);
   
end
%
rootfol = '../figures/';
saveas(ax,[rootfol 'Boxchart-NutsOnly-All-Glaciers-yearsgrouped.png']);

% 