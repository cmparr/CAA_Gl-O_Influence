function plotBox(data, lims, label, class)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot boxplots of 2 classes (glacerized vs unglacierized, etc)
%   class choices: glaciers, fjords, north
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

years = 2019:1:2022; alph = 0.10;

%groupsg = {'Glacierized';'Glacierized';'Glacierized';'Glacierized'; 'Unglacierized';'Unglacierized';'Unglacierized';'Unglacierized'};
groupsg = {'Glacierized';'Glacierized';'Glacierized';'Glacierized'; 'Riverine';'Riverine';'Riverine'};
%groupsg = {'Glacierized';'Riverine';'Riverine'};
groupsf = {'Open';'Open';'Open';'Fjorded'; 'Fjorded';'Fjorded';'Open' };
%groupsf = {'Open';'Open';'Open';'Fjorded'; 'Fjorded';'Fjorded';'Fjorded';'Open' };
groupsn = {'South';'South';'North';'North';'North'; 'North';'South'};

if strcmp(class, 'all')
groups = cat(2, groupsg, groupsf, groupsn);

% plotting
tiledlayout(4, 3, "TileSpacing","compact");

for ii = 1:4
    y = string(years(ii));
   %subplot(4,1,ii);
   nexttile
boxplot(reshape(data(ii,:,:),size(data, [2]),16)', groups(:, 1)); hold on;
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),'b','FaceAlpha',alph);
end
ylim(lims);
title(y); 
ylabel(label); 
%
nexttile
boxplot(reshape(data(ii,:,:),size(data, [2]),16)', groups(:, 2)); hold on;
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),'b','FaceAlpha',alph);
end
ylim(lims);
title(y);
%
nexttile
boxplot(reshape(data(ii,:,:),size(data, [2]),16)', groups(:, 3)); hold on;
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),'b','FaceAlpha',alph);
end
ylim(lims);
title(y);
end

else
    
if strcmp(class, 'glaciers')
    groups = groupsg;
elseif strcmp(class, 'fjords')
    groups = groupsf;
elseif strcmp(class, 'north')
    groups = groupsn;
else
    disp('unknown class')
    return
end

tiledlayout(4, 1, "TileSpacing","compact");

for ii = 1:4
    y = string(years(ii));
nexttile;
bh = boxplot(reshape(data(ii,:,:),size(data, [2]),16)', groups); hold on;
set(gca,'FontSize',13);
set(bh,'LineWidth', 1.1);
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),'k','FaceAlpha',alph);
end
ylim(lims);
title(y, 'Fontsize', 15); 
ylabel(label, 'FontSize',13); 
end

end
