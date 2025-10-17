%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot upwelling index with different variables (nitrate, chla, etc.)
% Figure 10. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear

set(0, 'defaultAxesFontSize', 13);
%
load('../Data/combined_nuts_ctd_all_noflag_dz5.mat');
% set g = 4 for only Jones Sound Glaciers 
glaciers = {'Jakeman','Sverdrup','Sydkap','Belcher','Trinity'};
gl_area = [498 765 491 1134 3085];
int_smb = zeros(length(glaciers),4); dz = int_smb; chla = dz; turb = dz; no3 = dz;
sample_day = [[210, nan, 240, nan];[216, 217, 225, 235]; [226, 227, 226, 234];[214, 212, 232, 222];[nan nan 215 nan]];

cols = cmocean('haline',4);
markers = {'pentagram';'o';'square';'^';'diamond'};

%% upwelling strength
for g = 1:4 %length(glaciers)
load('../Data/UpwellingIndex.mat');
T = T(T.dist <= 20,:);
T.Location = string(T.Location); T.Transect = string(T.Transect);
if g == 3
    T = T(contains([T.Location], {glaciers{g}, 'South'})&contains([T.Transect], 'out'),:);
    % for south cape fiord
elseif g == 5
    T = T(contains([T.Location], {'Talbot'})&contains([T.Transect], 'out'),:);
else
    T = T(contains([T.Location], {glaciers{g}})&contains([T.Transect], 'out'),:);
end
dz(g,:) = [mean(T(T.Year == 2019,:).idx, 'omitnan'); mean(T(T.Year == 2020,:).idx, 'omitnan'); mean(T(T.Year == 2021,:).idx, 'omitnan'); mean(T(T.Year == 2022,:).idx, 'omitnan')];
chla(g,:) = [mean(T(T.Year == 2019,:).chla, 'omitnan'); mean(T(T.Year == 2020,:).chla, 'omitnan'); mean(T(T.Year == 2021,:).chla, 'omitnan'); mean(T(T.Year == 2022,:).chla, 'omitnan')];
no3(g,:) = [mean(T(T.Year == 2019,:).no3, 'omitnan'); mean(T(T.Year == 2020,:).no3, 'omitnan'); mean(T(T.Year == 2021,:).no3, 'omitnan'); mean(T(T.Year == 2022,:).no3, 'omitnan')];
end

%% NO3 vs Chla
figure;
tiledlayout(1,2,"TileSpacing","tight","Padding","compact");
nexttile(2); hold on;
%yline(0, 'HandleVisibility','off');
%xline(0, 'HandleVisibility','off');
for g = 1:4 %length(glaciers)
    for yy = 1:4
        scatter(no3(g, yy), chla(g,yy), nh3(g,yy), 'filled', 'marker', markers{g}, ...
            'MarkerFaceColor',cols(yy,:), 'MarkerEdgeColor','k', 'HandleVisibility','off', 'MarkerFaceAlpha',0.85); hold on;
    end
end
xlabel('Mean Nitrate Concentration (\muM)'); xlim([-0.1 9]); 
ylabel('Mean Chlorophyll-\alpha Concentration (RFU)'); ylim([0 15]);
grid('minor');
for gg = 1:4 %5
    scatter(0, -10, 'filled','marker', markers{gg}, 'MarkerFaceColor',[1 1 1]*0.5,'MarkerEdgeColor','k');
end
for yy = 1:4
        scatter(0, -10, 'filled','marker', 'o', 'MarkerFaceColor',cols(yy,:),'MarkerEdgeColor','k');
end
%legend({'Jakeman';'Sverdrup';'Sydkap';'Belcher';'2019';'2020';'2021';'2022'}, 'NumColumns',2, 'Location','northwest');

%% idx vs no3
nexttile(1); hold on;
yline(0, 'HandleVisibility','off');
xline(0, 'HandleVisibility','off');
for g = 1: 4%length(glaciers)
    for yy = 1:4
        scatter(dz(g, yy), no3(g,yy), 80, 'filled', 'marker', markers{g}, ...
            'MarkerFaceColor',cols(yy,:), 'MarkerEdgeColor','k', 'HandleVisibility','off', 'MarkerFaceAlpha',0.85); hold on;
    end
end
xlabel('Upwelling Indicator'); xlim([-2 3.5]); 
ylabel('Mean Nitrate Concentration (\muM)'); ylim([-0.1 9]);
grid('minor');
for gg = 1:4 %length(glaciers)
    scatter(0, -10, 'filled','marker', markers{gg}, 'MarkerFaceColor',[1 1 1]*0.5,'MarkerEdgeColor','k');
end
for yy = 1:4
        scatter(0, -10, 'filled','marker', 'o', 'MarkerFaceColor',cols(yy,:),'MarkerEdgeColor','k');
end
legend({'Jakeman';'Sverdrup';'Sydkap';'Belcher';'2019';'2020';'2021';'2022'}, 'NumColumns',2, 'Location','northwest');
%legend({'Jakeman';'Sverdrup';'Sydkap';'Belcher';'2019';'2020';'2021';'2022'}, 'NumColumns',2, 'Location','northwest');

save(gca, '../figures/idxvsno3_no3vschla.png')