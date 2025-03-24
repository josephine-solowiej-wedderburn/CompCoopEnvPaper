filedr_savePlots = '/Users/josephinesolowiej-wedderburn/Documents/EnvProject/Plots/PiesForSingleRemoval/';

% %database
% collection = 'AGORA';
% oldCol = 'Agora';
% colAbr = 'AG';
% num_ec = 424;
collection = 'CarveMe';
oldCol = 'CarveMe';
colAbr = 'CM';
num_ec = 499;

npairs = 1000; %rather than looking at all possible pairs, choose a set number

%load data (from useKR100GE data)
filedrdat =  '/Users/josephinesolowiej-wedderburn/Documents/EnvProject/Data_PaperFigures/';
fnam_onlyU = [filedrdat, 'HPC2N_output/', collection, '_SingleRemoval_1to2500'];
%load pairs:
p_vec = 1:2500;
load([filedrdat, colAbr, '10000pairs.mat']);
%load summary stats to see no always essential
load([filedrdat, collection, '_bigRun100GEsCombined_summary_data_ECs100_conc1000.mat']);
%load the single removal data summarised across the pairs
filename_as_string = [fnam_onlyU, '/Allps', 'n_used_coop.mat'];
load(filename_as_string)
filename_as_string = [fnam_onlyU, '/Allps', 'n_used_comp.mat'];
load(filename_as_string)
%categories when removing one by one
filename_as_string = [fnam_onlyU, '/Allps', 'coop_newcatcount.mat'];
load(filename_as_string)
filename_as_string = [fnam_onlyU, '/Allps', 'comp_newcatcount.mat'];
load(filename_as_string)
%env tracker when removing one by one
filename_as_string = [fnam_onlyU, '/Allps', 'EnvTrackerCoop.mat'];
load(filename_as_string)
filename_as_string = [fnam_onlyU, '/Allps', 'EnvTrackerComp.mat'];
load(filename_as_string)
%did they comp/coop
filename_as_string = [fnam_onlyU, '/Allps', 'didTheyFCoop.mat'];
load(filename_as_string)
filename_as_string = [fnam_onlyU, '/Allps', 'didTheyComp.mat'];
load(filename_as_string)

%% 
%not all pairs competed and/or cooperated; restrict the search to the pairs
%that did both for a fair comparison?

candcpairs = p_vec(didTheyFCoop==1 & didTheyComp==1);

candcpairs = datasample(candcpairs, npairs, 'Replace', false);

coop_newcatcount = coop_newcatcount(candcpairs,:);
comp_newcatcount = comp_newcatcount(candcpairs,:);

%% count stats

n_alwaysEss = summary_data_ECs100_conc1000(candcpairs, 22);

cac_usedComp = n_used_comp(candcpairs) - n_alwaysEss;
cac_usedCoop = n_used_coop(candcpairs) - n_alwaysEss;

n_coopStaycoop = coop_newcatcount(:,1);
n_compStaycomp = comp_newcatcount(:,4);

n_coop2comp = coop_newcatcount(:,4);
n_comp2coop = comp_newcatcount(:,1);

n_coop21wO = (coop_newcatcount(:,10) + coop_newcatcount(:,13));
n_comp21wO = (comp_newcatcount(:,10) + comp_newcatcount(:,13));

n_coop22wO = coop_newcatcount(:,16);
n_comp22wO = comp_newcatcount(:,16);

n_coop2newEss = coop_newcatcount(:,17);
n_comp2newEss = comp_newcatcount(:,17);

n_coop2other_v1 = coop_newcatcount(:,2) + coop_newcatcount(:,3) + coop_newcatcount(:,5) + coop_newcatcount(:,6) + coop_newcatcount(:,7) + coop_newcatcount(:,8) + coop_newcatcount(:,9) + coop_newcatcount(:,11) + coop_newcatcount(:,12) + coop_newcatcount(:,14) + coop_newcatcount(:,15);
n_comp2other_v1 = comp_newcatcount(:,2) + comp_newcatcount(:,3) + comp_newcatcount(:,5) + comp_newcatcount(:,6) + comp_newcatcount(:,7) + comp_newcatcount(:,8) + comp_newcatcount(:,9) + comp_newcatcount(:,11) + comp_newcatcount(:,12) + comp_newcatcount(:,14) + comp_newcatcount(:,15);

% means:
mean_cac_usedCoop = mean(cac_usedCoop);
mean_cac_usedComp = mean(cac_usedComp);

mean_coopStaycoop = mean(n_coopStaycoop);
mean_compStaycomp = mean(n_compStaycomp);

mean_coop2comp = mean(n_coop2comp);
mean_comp2coop = mean(n_comp2coop);

mean_coop21wO = mean(n_coop21wO);
mean_comp21wO = mean(n_comp21wO);

mean_coop22wO = mean(n_coop22wO);
mean_comp22wO = mean(n_comp22wO);

mean_coop2newEss = mean(n_coop2newEss);
mean_comp2newEss = mean(n_comp2newEss);

mean_coop2other_v1 = mean(n_coop2other_v1);
mean_comp2other_v1 = mean(n_comp2other_v1);

% checkMeanCoop = (mean_coopStaycoop+mean_coop2comp+mean_coop21wO+mean_coop22wO+mean_coop2newEss+mean_coop2other_v1);
% checkMeanComp = (mean_compStaycomp+mean_comp2coop+mean_comp21wO+mean_comp22wO+mean_comp2newEss+mean_comp2other_v1);

%% test to see whether we have sampled enough:

%bin the data:
nbins = 10;
binsize = npairs/nbins;
whichPairs = datasample(1:npairs, nbins*binsize, 'Replace', false);
binSplit = repmat(1:nbins, 1, binsize);

% [p,tbl,stats] = anova1(cac_usedComp, binSplit);
% [p,tbl,stats] = anova1(cac_usedCoop, binSplit);
% [p,tbl,stats] = anova1(n_coopStaycoop, binSplit);
% [p,tbl,stats] = anova1(n_compStaycomp, binSplit);
% [p,tbl,stats] = anova1(n_coop2comp, binSplit);
% [p,tbl,stats] = anova1(n_comp2coop, binSplit);
% [p,tbl,stats] = anova1(n_coop21wO, binSplit);
% [p,tbl,stats] = anova1(n_comp21wO, binSplit);
% [p,tbl,stats] = anova1(n_coop22wO, binSplit);
% [p,tbl,stats] = anova1(n_comp22wO, binSplit);
% [p,tbl,stats] = anova1(n_coop2newEss, binSplit);
% [p,tbl,stats] = anova1(n_comp2newEss, binSplit);
% [p,tbl,stats] = anova1(n_coop2other_v1, binSplit);
% [p,tbl,stats] = anova1(n_comp2other_v1, binSplit);

%% permutation test to compare data from comp vs coop starting environments:

matcomp = [((100-cac_usedComp)+n_compStaycomp), n_comp2coop, (n_comp21wO + n_comp22wO), n_comp2other_v1, n_comp2newEss];
matcoop = [((100-cac_usedCoop)+n_coopStaycoop), n_coop2comp, (n_coop21wO + n_coop22wO), n_coop2other_v1, n_coop2newEss];

% Compute observed mean difference (Euclidean distance)
observed_diff = norm(mean(matcomp) - mean(matcoop));

% Permutation test setup
num_permutations = 10000;
perm_diffs = zeros(num_permutations, 1);

% Combine both datasets
combined = [matcomp; matcoop];
num_samples = size(combined, 1);

for i = 1:num_permutations
    shuffled = combined(randperm(num_samples), :); % Shuffle rows
    perm_diffs(i) = norm(mean(shuffled(1:size(matcomp,1), :)) - mean(shuffled(size(matcoop,1)+1:end, :)));
end

% Compute p-value
p = mean(perm_diffs >= observed_diff);
disp(['Permutation test p-value: ', num2str(p)]);

% 
% [p,tbl,stats] = anova1([((100-cac_usedComp)+n_compStaycomp);((100-cac_usedCoop)+n_coopStaycoop)], [ones(npairs,1);2*ones(npairs,1)]);
% [p,tbl,stats] = anova1([n_comp2coop;n_coop2comp], [ones(npairs,1);2*ones(npairs,1)]);
% [p,tbl,stats] = anova1([(n_comp21wO + n_comp22wO);(n_coop21wO + n_coop22wO)], [ones(npairs,1);2*ones(npairs,1)]);
% [p,tbl,stats] = anova1([n_comp2other_v1;n_coop2other_v1], [ones(npairs,1);2*ones(npairs,1)]);
% [p,tbl,stats] = anova1([n_comp2newEss;n_coop2newEss], [ones(npairs,1);2*ones(npairs,1)]);

%% define colours (as in fig 6)
col_fcoop = [184/247 225/247 134/247];
col_obl = [27/247 120/247 55/247 ];
col_oth = [0.75 0.75 0.75];    %maybe change to white or a lighter shade of grey?
col_comp = [208/247 28/247 139/247];
col_noGrow = [0 0 0];

%% make pies for the fraction ammounts:

%COMP STAY COMP
figure
ax = gca(); 
%no labels
labels = {''};
temp = floor((100-mean_cac_usedComp)+mean_compStaycomp);
frac = ((100-mean_cac_usedComp)+mean_compStaycomp) - temp;
p = pie(ax, frac, labels);
%change colours
cols = [col_comp ];
ax.Colormap = cols;
% %save
% figstr = [filedr_savePlots, colAbr, '_CompStayComp.fig'];
% saveas(gcf,figstr)
% % pngstr = [filedr_savePlots, colAbr, '_CompStayComp.png'];
% % saveas(gcf,pngstr)
% svgstr = [filedr_savePlots, colAbr, '_CompStayComp.svg'];
% saveas(gcf,svgstr)

%COMP 2 COOP
figure
ax = gca(); 
%no labels
labels = {''};
temp = floor(mean_comp2coop);
frac = mean_comp2coop - temp;
p = pie(ax, frac, labels);
%change colours
cols = [col_fcoop ];
ax.Colormap = cols;
% %save
% figstr = [filedr_savePlots, colAbr, '_Comp2Fcoop.fig'];
% saveas(gcf,figstr)
% % svgstr = [filedr_savePlots, colAbr, '_Comp2Fcoop.svg'];
% % saveas(gcf,svgstr)
% epsstr = [filedr_savePlots, colAbr, '_Comp2Fcoop.eps'];
% saveas(gcf,epsstr)

%COMP 2 OBL
figure
ax = gca(); 
%no labels
labels = {''};
temp = floor(mean_comp21wO + mean_comp22wO);
frac = (mean_comp21wO + mean_comp22wO) - temp;
p = pie(ax, frac, labels);
%change colours
cols = [col_obl ];
ax.Colormap = cols;
% %save
% figstr = [filedr_savePlots, colAbr, '_Comp2Obl.fig'];
% saveas(gcf,figstr)
% % svgstr = [filedr_savePlots, colAbr, '_Comp2Obl.svg'];
% % saveas(gcf,svgstr)
% epsstr = [filedr_savePlots, colAbr, '_Comp2Obl.eps'];
% saveas(gcf,epsstr)

%COMP 2 OTH
figure
ax = gca(); 
%no labels
labels = {''};
temp = floor(mean_comp2other_v1);
frac = mean_comp2other_v1 - temp;
p = pie(ax, frac, labels);
%change colours
cols = [col_oth ];
ax.Colormap = cols;
% %save
% figstr = [filedr_savePlots, colAbr, '_Comp2Oth.fig'];
% saveas(gcf,figstr)
% svgstr = [filedr_savePlots, colAbr, '_Comp2Oth.svg'];
% saveas(gcf,svgstr)

%COMP 2 NO GROW
figure
ax = gca(); 
%no labels
labels = {''};
temp = floor(mean_comp2newEss);
frac = mean_comp2newEss - temp;
p = pie(ax, frac, labels);
%change colours
cols = [col_noGrow ];
ax.Colormap = cols;
% %save
% figstr = [filedr_savePlots, colAbr, '_Comp2noGrow.fig'];
% saveas(gcf,figstr)
% svgstr = [filedr_savePlots, colAbr, '_Comp2noGrow.svg'];
% saveas(gcf,svgstr)

%COOP STAY COOP
figure
ax = gca(); 
%no labels
labels = {''};
temp = floor((100-mean_cac_usedCoop)+mean_coopStaycoop);
frac = ((100-mean_cac_usedCoop)+mean_coopStaycoop) - temp;
p = pie(ax, frac, labels);
%change colours
cols = [col_fcoop ];
ax.Colormap = cols;
% %save
% figstr = [filedr_savePlots, colAbr, '_CoopStayCoop.fig'];
% saveas(gcf,figstr)
% % pngstr = [filedr_savePlots, colAbr, '_CoopStayCoop.png'];
% % saveas(gcf,pngstr)
% svgstr = [filedr_savePlots, colAbr, '_CoopStayCoop.svg'];
% saveas(gcf,svgstr)

%COOP 2 COMP
figure
ax = gca(); 
%no labels
labels = {''};
temp = floor(mean_coop2comp);
frac = mean_coop2comp - temp;
p = pie(ax, frac, labels);
%change colours
cols = [col_comp ];
ax.Colormap = cols;
% %save
% figstr = [filedr_savePlots, colAbr, '_Coop2comp.fig'];
% saveas(gcf,figstr)
% svgstr = [filedr_savePlots, colAbr, '_Coop2comp.svg'];
% saveas(gcf,svgstr)

%COOP 2 OBL
figure
ax = gca(); 
%no labels
labels = {''};
temp = floor(mean_coop21wO + mean_coop22wO);
frac = (mean_coop21wO + mean_coop22wO) - temp;
p = pie(ax, frac, labels);
%change colours
cols = [col_obl ];
ax.Colormap = cols;
% %save
% figstr = [filedr_savePlots, colAbr, '_Coop2Obl.fig'];
% saveas(gcf,figstr)
% svgstr = [filedr_savePlots, colAbr, '_Coop2Obl.svg'];
% saveas(gcf,svgstr)

%COOP 2 OTH
figure
ax = gca(); 
%no labels
labels = {''};
temp = floor(mean_coop2other_v1);
frac = mean_coop2other_v1 - temp;
p = pie(ax, frac, labels);
%change colours
cols = [col_oth ];
ax.Colormap = cols;
% %save
% figstr = [filedr_savePlots, colAbr, '_Coop2Oth.fig'];
% saveas(gcf,figstr)
% svgstr = [filedr_savePlots, colAbr, '_Coop2Oth.svg'];
% saveas(gcf,svgstr)

%COMP 2 NO GROW
figure
ax = gca(); 
%no labels
labels = {''};
temp = floor(mean_coop2newEss);
frac = mean_coop2newEss - temp;
p = pie(ax, frac, labels);
%change colours
cols = [col_noGrow ];
ax.Colormap = cols;
% %save
% figstr = [filedr_savePlots, colAbr, '_Coop2noGrow.fig'];
% saveas(gcf,figstr)
% svgstr = [filedr_savePlots, colAbr, '_Coop2noGrow.svg'];
% saveas(gcf,svgstr)


%% histograms:

% %%start in f coop
% 
%histogram: facultative cooperation -> competition
figure
% h_fc2c = histogram(n_coop2comp, 'Normalization','percentage');
edges = [-0.25:0.5:22.25];
h_fc2c = histogram(n_coop2comp,edges, 'Normalization','percentage');
hold on
% h_fc2c.FaceColor = [208,28,139]/247;
h_fc2c.FaceColor = col_comp;
h_fc2c.BinWidth = 1;
h_fc2c.FaceAlpha  = 1;
xlim([-0.25,12])
ylim([0,60])
% xlabel('\# compounds switch fac. coop to comp','Interpreter','Latex')
% ylabel('% pairs')
set(gca, 'FontSize', 50)
set(gca, 'FontName', 'Times New Roman')
set(gca,'LineWidth',3,'TickLabelInterpreter','Latex','TickLength',[.01 .01],'TickDir','out','Box','off');
pbaspect([4 1 1])
%save
figstr = [filedr_savePlots, colAbr, '_SingRemHISTDiscrete_fc2c.fig'];
saveas(gcf,figstr)
pngstr = [filedr_savePlots, colAbr, '_SingRemHISTDiscrete_fc2c.png'];
saveas(gcf,pngstr)

%histogram: facultative cooperation -> obligate
figure
% h_fc2o = histogram(n_coop21wO+n_coop22wO, 'Normalization','percentage');
edges = [-0.25:0.5:22.25];
h_fc2o = histogram(n_coop21wO+n_coop22wO, edges, 'Normalization','percentage');
hold on
h_fc2o.FaceColor = col_obl;
h_fc2o.BinWidth = 1;
h_fc2o.FaceAlpha  = 1;
xlim([-0.25,12])
ylim([0,60])
% xlabel('\# compounds switch fac. coop to obligate','Interpreter','Latex')
% ylabel('% pairs')
set(gca, 'FontSize', 50)
set(gca, 'FontName', 'Times New Roman')
set(gca,'LineWidth',3,'TickLabelInterpreter','Latex','TickLength',[.01 .01],'TickDir','out','Box','off');
pbaspect([4 1 1])
%save
figstr = [filedr_savePlots, colAbr, '_SingRemHISTDiscrete_fc2o.fig'];
saveas(gcf,figstr)
pngstr = [filedr_savePlots, colAbr, '_SingRemHISTDiscrete_fc2o.png'];
saveas(gcf,pngstr)% %variable bin widths
 
% %%start in comp

%histogram: competition -> facultative cooperation
figure
% h_c2fc = histogram(n_comp2coop, 'Normalization','percentage');
edges = [-0.25:0.5:22.25];
h_c2fc = histogram(n_comp2coop, edges, 'Normalization','percentage');
hold on
h_c2fc.FaceColor = col_fcoop;
h_c2fc.BinWidth = 1;
h_c2fc.FaceAlpha  = 1;
xlim([-0.25,12])
ylim([0,60])
% xlabel('\# compounds switch comp to fac. coop','Interpreter','Latex')
% ylabel('% pairs')
set(gca, 'FontSize', 50)
set(gca, 'FontName', 'Times New Roman')
set(gca,'LineWidth',3,'TickLabelInterpreter','Latex','TickLength',[.01 .01],'TickDir','out','Box','off');
pbaspect([4 1 1])
%save
figstr = [filedr_savePlots, colAbr, '_SingRemHISTDiscrete_c2fc.fig'];
saveas(gcf,figstr)
pngstr = [filedr_savePlots, colAbr, '_SingRemHISTDiscrete_c2fc.png'];
saveas(gcf,pngstr)
 
%histogram: competition -> obligate
figure
% h_c2o = histogram(n_comp21wO+n_comp22wO, 'Normalization','percentage');
edges = [-0.25:0.5:22.25];
h_c2o = histogram(n_comp21wO+n_comp22wO, edges, 'Normalization','percentage');
hold on
h_c2o.FaceColor = col_obl;
h_c2o.BinWidth = 1;
h_c2o.FaceAlpha  = 1;
xlim([-0.25,12])
ylim([0,60])
% xlabel('\# compounds switch comp to obligate','Interpreter','Latex')
% ylabel('% pairs')
set(gca, 'FontSize', 50)
set(gca, 'FontName', 'Times New Roman')
set(gca,'LineWidth',3,'TickLabelInterpreter','Latex','TickLength',[.01 .01],'TickDir','out','Box','off');
pbaspect([4 1 1])
%save
figstr = [filedr_savePlots, colAbr, '_SingRemHISTDiscrete_c2o.fig'];
saveas(gcf,figstr)
pngstr = [filedr_savePlots, colAbr, '_SingRemHISTDiscrete_c2o.png'];
saveas(gcf,pngstr)


