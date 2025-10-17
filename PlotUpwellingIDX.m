%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot upwelling index with different variables (nitrate, chla, etc.)
% Figure 10. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear

set(0, 'defaultAxesFontSize', 13);
%
load('C:\Users\claire\Dropbox (Bhatia Lab)\Bhatia Lab Team Folder\Projects\Devon\Multi-year_analyses\Analysis\nutrients\master_sheets\combined_CTD\combined_nuts_ctd_all_noflag_dz5.mat');
% set g = 4 for only Jones Sound Glaciers 
glaciers = {'Jakeman','Sverdrup','Sydkap','Belcher','Trinity'};
gl_area = [498 765 491 1134 3085];
int_smb = zeros(length(glaciers),4); dz = int_smb; chla = dz; turb = dz; no3 = dz;
sample_day = [[210, nan, 240, nan];[216, 217, 225, 235]; [226, 227, 226, 234];[214, 212, 232, 222];[nan nan 215 nan]];

cols = cmocean('haline',4);
markers = {'pentagram';'o';'square';'^';'diamond'};

%% Annual integrated runoff
for g = 1:length(glaciers)
load(['C:\Users\claire\Dropbox (Bhatia Lab)\Bhatia Lab Team Folder\Projects\Devon\Multi-year_analyses\Mass_Balance\Bhatia\Output_matfiles\' glaciers{g} '-DailyRACMOOutput-2015-2022.mat']); 
smb = reshape(glacier_struct.smb',[],1); smb(smb >= 0) = 0;
for yy = 2019:2022
    %{
    % total annual runoff before sampling day
    smby = smb(glacier_struct.time > datetime(yy,1,1) & glacier_struct.time < datetime(yy,1,sample_day(g,yy-2018)));
    smby = -smby;
    int_smb(g, yy-2018) = sum(smby, 'omitnan')*sample_day(g,yy-2018)*gl_area(g)/1e12; %%% Following Williams 2021 takes negative SMB * Gl.Area * N days to end of year / 1e12 (convert cm to m, I guess smb is in cm?) 
    %}
    % runoff on sampling day 
    sample = datetime(yy,1,sample_day(g,yy-2018));
    if isnat(sample)
    int_smb(g,yy-2018) = NaN;
    else
    int_smb(g,yy-2018) = mean(-smb(glacier_struct.time >= sample - 1 & glacier_struct.time <= sample + 1)*gl_area(g))/1e6;
    end
    %}
end
end

%% upwelling strength
for g = 1:4 %length(glaciers)
load('C:\Users\claire\Dropbox (Bhatia Lab)\Bhatia Lab Team Folder\Projects\Devon\Multi-year_analyses\Analysis\nutrients\master_sheets\combined_CTD\UpwellingIndex.mat');
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
turb(g,:) = [mean(T(T.Year == 2019,:).turb, 'omitnan'); mean(T(T.Year == 2020,:).turb, 'omitnan'); mean(T(T.Year == 2021,:).turb, 'omitnan'); mean(T(T.Year == 2022,:).turb, 'omitnan')];
nh3(g,:) = [mean(T(T.Year == 2019,:).NH3, 'omitnan'); mean(T(T.Year == 2020,:).NH3, 'omitnan'); mean(T(T.Year == 2021,:).NH3, 'omitnan'); mean(T(T.Year == 2022,:).NH3, 'omitnan')];
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
%xlim([-1 max(int_smb,[],'all')+0.1*max(int_smb,[],'all')]);
%xlabel('Total Discharge Volume (Gt/Yr)'); xlim([0 max(int_smb,[],'all')+0.1*max(int_smb,[],'all')]);
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
%xlim([-1 max(int_smb,[],'all')+0.1*max(int_smb,[],'all')]);
%xlabel('Total Discharge Volume (Gt/Yr)'); xlim([0 max(int_smb,[],'all')+0.1*max(int_smb,[],'all')]);
grid('minor');
for gg = 1:4 %length(glaciers)
    scatter(0, -10, 'filled','marker', markers{gg}, 'MarkerFaceColor',[1 1 1]*0.5,'MarkerEdgeColor','k');
end
for yy = 1:4
        scatter(0, -10, 'filled','marker', 'o', 'MarkerFaceColor',cols(yy,:),'MarkerEdgeColor','k');
end
legend({'Jakeman';'Sverdrup';'Sydkap';'Belcher';'2019';'2020';'2021';'2022'}, 'NumColumns',2, 'Location','northwest');
%legend({'Jakeman';'Sverdrup';'Sydkap';'Belcher';'Talbot';'2019';'2020';'2021';'2022'}, 'NumColumns',2, 'Location','northwest');

%% 3D scatter
%
figure;
for g = 1:length(glaciers)
    for yy = 1:4
        scatter3(turb(g, yy), no3(g,yy), chla(g,yy), 80, 'filled', 'marker', markers{g}, ...
            'MarkerFaceColor',cols(yy,:), 'MarkerEdgeColor','k', 'HandleVisibility','off'); hold on;
    end
end
%
xlabel('Mean Turbidity upper 30 m (NTU)'); xlim([0 15]); 
%xlabel('Upwelling Index'); xlim([-6.5 3.5]); 
ylabel('Mean nitrate concentration (\muM)'); ylim([0 9]);
zlabel('Mean chlorophyll-\alpha concentration (\mug/L)'); zlim([0 15]);
%xlim([-1 max(int_smb,[],'all')+0.1*max(int_smb,[],'all')]);
%xlabel('Total Discharge Volume (Gt/Yr)'); xlim([0 max(int_smb,[],'all')+0.1*max(int_smb,[],'all')]);
grid('minor');
for gg = 1:length(glaciers)
    scatter(0, -10, 'filled','marker', markers{gg}, 'MarkerFaceColor',[1 1 1]*0.5,'MarkerEdgeColor','k');
end
for yy = 1:4
        scatter(0, -10, 'filled','marker', 'o', 'MarkerFaceColor',cols(yy,:),'MarkerEdgeColor','k');
end
legend({'Jakeman';'Sverdrup';'Sydkap';'Belcher';'2019';'2020';'2021';'2022'}, 'NumColumns',2, 'Location','northwest');
%}

%% Gl. dischg vs idx
figure; hold on;
yline(0, 'HandleVisibility','off');
xline(0, 'HandleVisibility','off');
for g = 1:length(glaciers)
    for yy = 1:4
        scatter(int_smb(g, yy), dz(g,yy), 80, 'filled', 'marker', markers{g}, ...
            'MarkerFaceColor',cols(yy,:), 'MarkerEdgeColor','k', 'HandleVisibility','off', 'MarkerFaceAlpha',0.85); hold on;
    end
end
ylabel('Upwelling Index'); ylim([-2 3.5]); 
%xlim([-1 max(int_smb,[],'all')+0.1*max(int_smb,[],'all')]);
xlabel('Sample Day Discharge (Kt/Yr)'); xlim([0 max(int_smb,[],'all')+0.1*max(int_smb,[],'all')]);
grid('minor');
for gg = 1:length(glaciers)
    scatter(0, -10, 'filled','marker', markers{gg}, 'MarkerFaceColor',[1 1 1]*0.5,'MarkerEdgeColor','k');
end
for yy = 1:4
        scatter(0, -10, 'filled','marker', 'o', 'MarkerFaceColor',cols(yy,:),'MarkerEdgeColor','k');
end
legend({'Jakeman';'Sverdrup';'Sydkap';'Belcher';'Talbot';'2019';'2020';'2021';'2022'}, 'NumColumns',2, 'Location','northwest');
