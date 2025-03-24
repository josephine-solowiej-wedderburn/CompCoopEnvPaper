%code to upload environments and measure the distance between comp envs,
%coop envs and all envs...

%% set up:

%collection:
% colAbr = 'AG';
colAbr = 'CM';

%processed data for FriendVsFoe project; upload from there
filedr = '/Users/josephinesolowiej-wedderburn/Documents/FriendVsFoeProject/data/';

%what size environments do we want to look at 
% nEC = 50;
nEC = 100;

%too many envs to measure everything, try:
nMeas = 1000;   %note crashes if trying to plot >~5000
n_bins = 20;

%% load env matrices:
%comp envs:
load([filedr, 'CompCoopEnvs_fromBigRun100GE/', colAbr, '_ECnum', num2str(nEC), 'conc1000_comp_envs.mat' ])

%f coop envs:
load([filedr, 'CompCoopEnvs_fromBigRun100GE/', colAbr, '_ECnum', num2str(nEC), 'conc1000_fcoop_envs.mat' ])

%obl envs:
load([filedr, 'OblEnvs_fromBigRun100GE/', colAbr, '_ECnum', num2str(nEC), 'conc1000_plusx_envs.mat' ])
load([filedr, 'OblEnvs_fromBigRun100GE/', colAbr, '_ECnum', num2str(nEC), 'conc1000_xplus_envs.mat' ])
load([filedr, 'OblEnvs_fromBigRun100GE/', colAbr, '_ECnum', num2str(nEC), 'conc1000_xx_envs.mat' ])


% %% calculate Jaccard and Hamming Distance
% 
% HD_compEnvs = zeros(n_bins,nchoosek(nMeas,2));
% HD_fcoopEnvs = zeros(n_bins,nchoosek(nMeas,2));
% HD_oblEnvs = zeros(n_bins,nchoosek(nMeas,2));
% HD_alllEnvs = zeros(n_bins,nchoosek(3*nMeas,2));
% 
% JD_compEnvs = zeros(n_bins,nchoosek(nMeas,2));
% JD_fcoopEnvs = zeros(n_bins,nchoosek(nMeas,2));
% JD_oblEnvs = zeros(n_bins,nchoosek(nMeas,2));
% JD_alllEnvs = zeros(n_bins,nchoosek(3*nMeas,2));
% 
% for ii = 1:n_bins
% 
%     disp(ii)
% 
%     whichRows = datasample(1:size(comp_envs,1), nMeas, 'Replace', false);
%     comp_envs_red = comp_envs(whichRows, :);
%     HD_compEnvs(ii,:) = pdist(comp_envs_red, 'hamming');
%     JD_compEnvs(ii,:) = pdist(comp_envs_red, 'jaccard');
% 
%     whichRows = datasample(1:size(fcoop_envs,1), nMeas, 'Replace', false);
%     fcoop_envs_red = fcoop_envs(whichRows, :);
%     HD_fcoopEnvs(ii,:) = pdist(fcoop_envs_red, 'hamming');
%     JD_fcoopEnvs(ii,:) = pdist(fcoop_envs_red, 'jaccard');
% 
%     whichRows = datasample(1:(size(plusx_envs,1)+size(xplus_envs,1)+size(xx_envs,1)), nMeas, 'Replace', false);
%     plusxRows = whichRows(whichRows<=size(plusx_envs,1));
%     xplusRows = whichRows(whichRows>size(plusx_envs,1) & whichRows<=(size(plusx_envs,1)+size(xplus_envs,1))) - size(plusx_envs,1);
%     xxRows = whichRows(whichRows>(size(plusx_envs,1)+size(xplus_envs,1))) - (size(plusx_envs,1)+size(xplus_envs,1));
%     obl_envs_red = [plusx_envs(plusxRows,:);
%                     xplus_envs(xplusRows,:);
%                     xx_envs(xxRows,:)];
%     HD_oblEnvs(ii,:) = pdist(obl_envs_red, 'hamming');
%     JD_oblEnvs(ii,:) = pdist(obl_envs_red, 'jaccard');
% 
%     HD_alllEnvs(ii,:) = pdist([comp_envs_red;fcoop_envs_red;obl_envs_red], 'hamming');
%     JD_alllEnvs(ii,:) = pdist([comp_envs_red;fcoop_envs_red;obl_envs_red], 'jaccard');
% 
% end
% 
% meanHD_compEnvs = mean(HD_compEnvs');
% meanHD_fcoopEnvs = mean(HD_fcoopEnvs');
% meanHD_oblEnvs = mean(HD_oblEnvs');
% meanHD_alllEnvs = mean(HD_alllEnvs');
% 
% meanJD_compEnvs = mean(JD_compEnvs');
% meanJD_fcoopEnvs = mean(JD_fcoopEnvs');
% meanJD_oblEnvs = mean(JD_oblEnvs');
% meanJD_alllEnvs = mean(JD_alllEnvs');
% 
% % %% plot HD:
% % figure
% % xHD = [ones(length(meanHD_compEnvs),1); 2*ones(length(meanHD_fcoopEnvs),1); 3*ones(length(meanHD_oblEnvs),1); 4*ones(length(meanHD_alllEnvs),1)];
% % yHD = [meanHD_compEnvs, meanHD_fcoopEnvs, meanHD_oblEnvs, meanHD_alllEnvs];
% % boxchart(xHD, yHD)
% % % 
% % %check normality:
% % normplot(meanHD_compEnvs)
% % normplot(meanHD_fcoopEnvs)
% % normplot(meanHD_oblEnvs)
% % normplot(meanHD_alllEnvs)
% % %
% % [pHD,tblHD,statsHD] = anova1(yHD, xHD);
% % figure
% % resultsHD = multcompare(statsHD);
% 
% %% plot JD:
% figure
% xJD = [ones(length(meanJD_compEnvs),1); 2*ones(length(meanJD_fcoopEnvs),1); 3*ones(length(meanJD_oblEnvs),1); 4*ones(length(meanJD_alllEnvs),1)];
% yJD = [meanJD_compEnvs, meanJD_fcoopEnvs, meanJD_oblEnvs, meanJD_alllEnvs];
% boxchart(xJD, yJD)
% % 
% %check normality:
% normplot(meanJD_compEnvs)
% normplot(meanJD_fcoopEnvs)
% normplot(meanJD_oblEnvs)
% normplot(meanJD_alllEnvs)
% %
% [pJD,tblJD,statsJD] = anova1(yJD, xJD);
% figure
% resultsJD = multcompare(statsJD);

%% calculate Jaccard Hamming Distance; lump together all coop envs

HD_compEnvs = zeros(n_bins,nchoosek(nMeas,2));
HD_coopEnvs = zeros(n_bins,nchoosek(nMeas,2));
HD_allEnvs = zeros(n_bins,nchoosek(2*nMeas,2));
HD_compvscoopEnvs = zeros(n_bins,nchoosek(nMeas,2));

JD_compEnvs = zeros(n_bins,nchoosek(nMeas,2));
JD_coopEnvs = zeros(n_bins,nchoosek(nMeas,2));
JD_allEnvs = zeros(n_bins,nchoosek(2*nMeas,2));
JD_compvscoopEnvs = zeros(n_bins,nchoosek(nMeas,2));

for ii = 1:n_bins
    
    disp(ii)

    whichRows = datasample(1:size(comp_envs,1), nMeas, 'Replace', false);
    comp_envs_red = comp_envs(whichRows, :);
    HD_compEnvs(ii,:) = pdist(comp_envs_red, 'hamming');
    JD_compEnvs(ii,:) = pdist(comp_envs_red, 'jaccard');

    whichRows = datasample(1:(size(plusx_envs,1)+size(xplus_envs,1)+size(xx_envs,1)+size(fcoop_envs,1)), nMeas, 'Replace', false);
    plusxRows = whichRows(whichRows<=size(plusx_envs,1));
    xplusRows = whichRows(whichRows>size(plusx_envs,1) & whichRows<=(size(plusx_envs,1)+size(xplus_envs,1))) - size(plusx_envs,1);
    xxRows = whichRows(whichRows>(size(plusx_envs,1)+size(xplus_envs,1)) & whichRows<=(size(plusx_envs,1)+size(xplus_envs,1)+size(xx_envs,1))) - (size(plusx_envs,1)+size(xplus_envs,1));
    fcoopoRows = whichRows(whichRows>(size(plusx_envs,1)+size(xplus_envs,1)+size(xx_envs,1))) - (size(plusx_envs,1)+size(xplus_envs,1)+size(xx_envs,1));
    coop_envs_red = [plusx_envs(plusxRows,:);
                    xplus_envs(xplusRows,:);
                    xx_envs(xxRows,:);
                    fcoop_envs(fcoopoRows, :)];
    HD_coopEnvs(ii,:) = pdist(coop_envs_red, 'hamming');
    JD_coopEnvs(ii,:) = pdist(coop_envs_red, 'jaccard');

    HD_allEnvs(ii,:) = pdist([comp_envs_red;coop_envs_red], 'hamming');
    JD_allEnvs(ii,:) = pdist([comp_envs_red;coop_envs_red], 'jaccard');

    for jj = 1:nMeas
        Onecomp_env = comp_envs_red(jj,:);
        HDsOnecomp_env = zeros(1,nMeas);
        JDsOnecomp_env = zeros(1,nMeas);
        for kk = 1:nMeas
            HDsOnecomp_env(kk) = pdist([Onecomp_env; coop_envs_red(kk,:)], 'hamming');
            JDsOnecomp_env(kk) = pdist([Onecomp_env; coop_envs_red(kk,:)], 'jaccard');
        end
        HD_compvscoopEnvs(ii, (nMeas*(jj-1)+1:nMeas*(jj))) = HDsOnecomp_env;
        JD_compvscoopEnvs(ii, (nMeas*(jj-1)+1:nMeas*(jj))) = JDsOnecomp_env;
    end

end

%% save data
filedr_save = '/Users/josephinesolowiej-wedderburn/Documents/EnvProject/Data_PaperFigures/SI/EnvSimilarity/';

fnam = [filedr_save, colAbr, 'EnvSize', num2str(nEC), '_', 'HD_compEnvs.mat'];
save(fnam, 'HD_compEnvs');
fnam = [filedr_save, colAbr, 'EnvSize', num2str(nEC), '_', 'HD_coopEnvs.mat'];
save(fnam, 'HD_coopEnvs');
fnam = [filedr_save, colAbr, 'EnvSize', num2str(nEC), '_', 'HD_allEnvs.mat'];
save(fnam, 'HD_allEnvs');
fnam = [filedr_save, colAbr, 'EnvSize', num2str(nEC), '_', 'HD_compvscoopEnvs.mat'];
save(fnam, 'HD_compvscoopEnvs');

fnam = [filedr_save, colAbr, 'EnvSize', num2str(nEC), '_', 'JD_compEnvs.mat'];
save(fnam, 'JD_compEnvs');
fnam = [filedr_save, colAbr, 'EnvSize', num2str(nEC), '_', 'JD_coopEnvs.mat'];
save(fnam, 'JD_coopEnvs');
fnam = [filedr_save, colAbr, 'EnvSize', num2str(nEC), '_', 'JD_allEnvs.mat'];
save(fnam, 'JD_allEnvs');
fnam = [filedr_save, colAbr, 'EnvSize', num2str(nEC), '_', 'JD_compvscoopEnvs.mat'];
save(fnam, 'JD_compvscoopEnvs');

%% calculate means
meanHD_compEnvs = mean(HD_compEnvs');
meanHD_coopEnvs = mean(HD_coopEnvs');
meanHD_alllEnvs = mean(HD_allEnvs');
meanHD_compvscoopEnvs = mean(HD_compvscoopEnvs');

meanJD_compEnvs = mean(JD_compEnvs');
meanJD_coopEnvs = mean(JD_coopEnvs');
meanJD_alllEnvs = mean(JD_allEnvs');
meanJD_compvscoopEnvs = mean(JD_compvscoopEnvs');

% %% plot HD:
% figure
% xHD = [ones(length(meanHD_compEnvs),1); 2*ones(length(meanHD_fcoopEnvs),1); 3*ones(length(meanHD_oblEnvs),1); 4*ones(length(meanHD_alllEnvs),1)];
% yHD = [meanHD_compEnvs, meanHD_fcoopEnvs, meanHD_oblEnvs, meanHD_alllEnvs];
% boxchart(xHD, yHD)
% % 
% %check normality:
% normplot(meanHD_compEnvs)
% normplot(meanHD_fcoopEnvs)
% normplot(meanHD_oblEnvs)
% normplot(meanHD_alllEnvs)
% %
% [pHD,tblHD,statsHD] = anova1(yHD, xHD);
% figure
% resultsHD = multcompare(statsHD);

%% plot JD:
%plot mean distance:
figure
xJD = [ones(length(meanJD_compEnvs),1); 2*ones(length(meanJD_coopEnvs),1); 3*ones(length(meanJD_compvscoopEnvs),1)];
name_order = convertCharsToStrings([{'comp'},{'coop'},{'comp vs coop'}]);
xJD_named = categorical(xJD, 1:3, name_order);
yJD = [meanJD_compEnvs, meanJD_coopEnvs, meanJD_compvscoopEnvs];
boxchart(xJD_named, yJD)
ylab = ylabel('mean Jaccard distance');
ylab.FontName = 'Times';
set(gca,'LineWidth',3,'TickLabelInterpreter','Latex','TickLength',[.015 .015],'FontSize',28,'Box','off');
axis square
%save
filedr_savePlots = '/Users/josephinesolowiej-wedderburn/Documents/EnvProject/Plots/SI_envSimilarity/';
figstr = [filedr_savePlots, colAbr, '_meanJD', num2str(n_bins), 'binsCompCoopVs_EnvSize', num2str(nEC), '.fig'];
saveas(gcf,figstr)
svgstr = [filedr_savePlots, colAbr, '_meanJD', num2str(n_bins), 'binsCompCoopVs_EnvSize', num2str(nEC), '.svg'];
saveas(gcf,svgstr)
pngstr = [filedr_savePlots, colAbr, '_meanJD', num2str(n_bins), 'binsCompCoopVs_EnvSize', num2str(nEC), '.png'];
saveas(gcf,pngstr)

% % %plot all the distances:
% % xallJD = repmat(1:n_bins,nchoosek(nMeas,2));
% % 
% %check normality:
% normplot(meanJD_compEnvs)
% normplot(meanJD_coopEnvs)
% normplot(meanJD_compvscoopEnvs)
% %
% [pJD,tblJD,statsJD] = anova1(yJD, xJD);
% figure
% resultsJD = multcompare(statsJD);
% 
% 
% 

%% load data and plot as bar graph instaad:

%%collection:
colAbr = 'AG';
% colAbr = 'CM';
%%what size environments do we want to look at 
% nEC = 50;
nEC = 100;
%%bins and size:
nMeas = 1000;
n_bins = 20;

filedr_save = '/Users/josephinesolowiej-wedderburn/Documents/EnvProject/Data_PaperFigures/SI/EnvSimilarity/';

load([filedr_save, colAbr, 'EnvSize', num2str(nEC), '_', 'JD_compEnvs.mat']);
load([filedr_save, colAbr, 'EnvSize', num2str(nEC), '_', 'JD_coopEnvs.mat']);
load([filedr_save, colAbr, 'EnvSize', num2str(nEC), '_', 'JD_compvscoopEnvs.mat']);

%bar graph:
figure
x = [1,2,3];
data = [ mean(mean(JD_compEnvs')), mean(mean(JD_coopEnvs')), mean(mean(JD_compvscoopEnvs'))]';
errhigh = [ std(mean(JD_compEnvs')./nMeas), std(mean(JD_coopEnvs')./nMeas), std(mean(JD_compvscoopEnvs')./nMeas)]';
b = bar(x,data,'FaceColor','flat'); 
b.FaceColor = [198,219,239]/255;
b.BarWidth = 0.6;
set(gca,'XTickLabel',{'comp','coop','comp vs coop'});
hold on
er = errorbar(x, data, errhigh, '.');
er.LineWidth = 3;
er.Color = [0 0 0];
axis square
set(gca,'LineWidth',3,'TickLabelInterpreter','Latex','TickLength',[.015 .015],'FontSize',28,'Box','off');
ylabel('mean Jaccard distance','FontSize',28,'Interpreter','Latex');
ylim([0,1])
xlim([0.5, 3.5])
h=gca; 
h.XAxis.TickLength = [0 0];

% %save
% filedr_savePlots = '/Users/josephinesolowiej-wedderburn/Documents/EnvProject/Plots/SI_envSimilarity/';
% figstr = [filedr_savePlots, colAbr, '_Bars_meanJD', num2str(n_bins), 'binsCompCoopVs_EnvSize', num2str(nEC), '.fig'];
% saveas(gcf,figstr)
% % svgstr = [filedr_savePlots, colAbr, '_Bars_meanJD', num2str(n_bins), 'binsCompCoopVs_EnvSize', num2str(nEC), '.svg'];
% % saveas(gcf,svgstr)
% pngstr = [filedr_savePlots, colAbr, '_Bars_meanJD', num2str(n_bins), 'binsCompCoopVs_EnvSize', num2str(nEC), '.png'];
% saveas(gcf,pngstr)
