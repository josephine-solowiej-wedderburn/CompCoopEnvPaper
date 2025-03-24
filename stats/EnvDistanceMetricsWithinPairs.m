%code to upload environments and measure the distance between comp envs,
%coop envs and all envs...

%NOW: only look at distance between environments of a pair, specify pairs
%that have at least 10 comp and 10 coop envs

%% set up:

%collection:
% collection = 'AGORA';
% colAbr = 'AG';
collection = 'CarveMe';
colAbr = 'CM';

% %processed data for FriendVsFoe project; upload from there
% filedr_extra = '/Users/josephinesolowiej-wedderburn/Documents/FriendVsFoeProject/data/';
%and data file for orig project:
filedr_data = '/Users/josephinesolowiej-wedderburn/Documents/EnvProject/Data_PaperFigures/';
filedr_dataBigRun = '/Users/josephinesolowiej-wedderburn/Documents/General_MetModels/Data/EnvProject_HPC2N_output/';
%path:
%use functions from 'General_MetModels' folder to get the desired envs
pathWcommonfuns = '/Users/josephinesolowiej-wedderburn/Documents/General_MetModels/Code/common_functions/';
addpath(pathWcommonfuns)

%what size environments do we want to look at 
nEC = 50;
% nEC = 100;

% %too many envs to measure everything, try:
% npairs = 1000; 

%min number of each type of environment that means we look at that pair:
min_envs = 5;

% %% load env matrices and pairs:
% %comp envs:
% load([filedr_extra, 'CompCoopEnvs_fromBigRun100GE/', colAbr, '_ECnum', num2str(nEC), 'conc1000_comp_envs.mat' ])
% load([filedr_extra, 'CompCoopEnvs_fromBigRun100GE/', colAbr, '_ECnum', num2str(nEC), 'conc1000_MiMj_comp_mat.mat' ])
% 
% %f coop envs:
% load([filedr_extra, 'CompCoopEnvs_fromBigRun100GE/', colAbr, '_ECnum', num2str(nEC), 'conc1000_fcoop_envs.mat' ])
% load([filedr_extra, 'CompCoopEnvs_fromBigRun100GE/', colAbr, '_ECnum', num2str(nEC), 'conc1000_MiMj_fcoop_mat.mat' ])
% 
% %obl envs:
% load([filedr_extra, 'OblEnvs_fromBigRun100GE/', colAbr, '_ECnum', num2str(nEC), 'conc1000_plusx_envs.mat' ])
% load([filedr_extra, 'OblEnvs_fromBigRun100GE/', colAbr, '_ECnum', num2str(nEC), 'conc1000_xplus_envs.mat' ])
% load([filedr_extra, 'OblEnvs_fromBigRun100GE/', colAbr, '_ECnum', num2str(nEC), 'conc1000_xx_envs.mat' ])
% load([filedr_extra, 'CompCoopEnvs_fromBigRun100GE/', colAbr, '_ECnum', num2str(nEC), 'conc1000_MiMj_plusx_mat.mat' ])
% load([filedr_extra, 'CompCoopEnvs_fromBigRun100GE/', colAbr, '_ECnum', num2str(nEC), 'conc1000_MiMj_xplus_mat.mat' ])
% load([filedr_extra, 'CompCoopEnvs_fromBigRun100GE/', colAbr, '_ECnum', num2str(nEC), 'conc1000_MiMj_xx_mat.mat' ])



%% load summary data to check which pairs have the right number of envs:
load([ filedr_data, colAbr, 'summary_data_ECs', num2str(nEC), '_conc1000_100GEsx4envs.mat' ])
eval(['summary_data = summary_data_ECs', num2str(nEC), '_conc1000;'])

%%JUST LOOK AT COMP AND F COOP:
types = ['(-/-)'; '(+/+)'];
%select pairs with enough
Allpair_vec = 1:size(summary_data, 1);
pairs_withEnough = Allpair_vec(summary_data(:,1)>=min_envs & summary_data(:,4)>=min_envs);

%randomly select npairs
% selectpairs_withEnough = datasample(pairs_withEnough, npairs, 'replace', false);
selectpairs_withEnough = pairs_withEnough;
npairs = length(selectpairs_withEnough); 

%% run through each selected pair, upload environments and calculate the hamming distances
%set up to store data:
% %environment hamming/jaccard distances:
% HD_compenvs = cell(npairs,1);
% HD_fcoopenvs = cell(npairs,1);
% HD_ALLenvs = cell(npairs,1);
% JD_compenvs = cell(npairs,1);
% JD_fcoopenvs = cell(npairs,1);
% JD_ALLenvs = cell(npairs,1);
%median and quartiles:
HDmedQ_compenvs = zeros(npairs,3);
HDmedQ_fcoopenvs = zeros(npairs,3);
HDmedQ_compvsfcoopenvs = zeros(npairs,3);
HDmedQ_ALLenvs = zeros(npairs,3);
JDmedQ_compenvs = zeros(npairs,3);
JDmedQ_fcoopenvs = zeros(npairs,3);
JDmedQ_compvsfcoopenvs = zeros(npairs,3);
JDmedQ_ALLenvs = zeros(npairs,3);

%mean:
HDmeanQ_compenvs = zeros(npairs,1);
HDmeanQ_fcoopenvs = zeros(npairs,1);
HDmeanQ_compvsfcoopenvs = zeros(npairs,1);
HDmeanQ_ALLenvs = zeros(npairs,1);
JDmeanQ_compenvs = zeros(npairs,1);
JDmeanQ_fcoopenvs = zeros(npairs,1);
JDmeanQ_compvsfcoopenvs = zeros(npairs,1);
JDmeanQ_ALLenvs = zeros(npairs,1);

for pp = 1:npairs

    output_envs = {};
    for ii = 1:size(types,1)
        int_type = types(ii,:);
        output_envs{ii} = fun_getEnvsFromBigRun100GEdata( selectpairs_withEnough(pp), collection, int_type, nEC, 1000 );
    end

    %comp envs
    temp_comp_envs = output_envs{1};
    temp_comp_envs = temp_comp_envs(:,3:end);
    temp_HD_compenvs = pdist(temp_comp_envs, 'hamming');
    temp_JD_compenvs = pdist(temp_comp_envs, 'jaccard');
    % HD_compenvs{pp} = temp_HD_compenvs;
    % JD_compenvs{pp} = temp_JD_compenvs;
    HDmedQ_compenvs(pp,:) = quantile(temp_HD_compenvs,[0.25 0.5 0.75]);
    JDmedQ_compenvs(pp,:) = quantile(temp_JD_compenvs,[0.25 0.5 0.75]);
    HDmeanQ_compenvs(pp) = mean(temp_HD_compenvs);
    JDmeanQ_compenvs(pp) = mean(temp_JD_compenvs);

    %f coop envs
    temp_fcoop_envs = output_envs{2};
    temp_fcoop_envs = temp_fcoop_envs(:,3:end);
    temp_HD_fcoopenvs = pdist(temp_fcoop_envs, 'hamming');
    temp_JD_fcoopenvs = pdist(temp_fcoop_envs, 'jaccard');
    % HD_compenvs{pp} = temp_HD_compenvs;
    % JD_compenvs{pp} = temp_JD_compenvs;
    HDmedQ_fcoopenvs(pp,:) = quantile(temp_HD_fcoopenvs,[0.25 0.5 0.75]);
    JDmedQ_fcoopenvs(pp,:) = quantile(temp_JD_fcoopenvs,[0.25 0.5 0.75]);
    HDmeanQ_fcoopenvs(pp) = mean(temp_HD_fcoopenvs); 
    JDmeanQ_fcoopenvs(pp) = mean(temp_JD_fcoopenvs);

    %compare each comp environment with all the f coop envs:
    n_comp = size(temp_comp_envs,1);
    n_fcoop = size(temp_fcoop_envs,1);
    temp_HDs_compvsfcoopenvs = zeros(n_comp,n_fcoop);
    temp_JDs_compvsfcoopenvs = zeros(n_comp,n_fcoop);
    for cc = 1:n_comp
        comp_vec = temp_comp_envs(cc,:);
        for ff = 1:n_fcoop
            temp_HDs_compvsfcoopenvs(cc,ff) = pdist([comp_vec; temp_fcoop_envs(ff,:)] , 'hamming');
            temp_JDs_compvsfcoopenvs(cc,ff) = pdist([comp_vec; temp_fcoop_envs(ff,:)] , 'jaccard');
        end
    end
    temp_HDs_compvsfcoopenvs = reshape(temp_HDs_compvsfcoopenvs, 1, []);
    HDmedQ_compvsfcoopenvs(pp,:) = quantile(temp_HDs_compvsfcoopenvs,[0.25 0.5 0.75]);
    HDmeanQ_compvsfcoopenvs(pp) = mean(temp_HDs_compvsfcoopenvs);
    temp_JDs_compvsfcoopenvs = reshape(temp_JDs_compvsfcoopenvs, 1, []);
    JDmedQ_compvsfcoopenvs(pp,:) = quantile(temp_JDs_compvsfcoopenvs,[0.25 0.5 0.75]);
    JDmeanQ_compvsfcoopenvs(pp) = mean(temp_JDs_compvsfcoopenvs);

    %all envs
    temp_all_envs = [temp_comp_envs; temp_fcoop_envs];
    temp_HD_allenvs = pdist(temp_all_envs, 'hamming');
    temp_JD_allenvs = pdist(temp_all_envs, 'jaccard');
    % HD_compenvs{pp} = temp_HD_compenvs;
    % JD_compenvs{pp} = temp_JD_compenvs;
    HDmedQ_ALLenvs(pp,:) = quantile(temp_HD_allenvs,[0.25 0.5 0.75]);
    JDmedQ_ALLenvs(pp,:) = quantile(temp_JD_allenvs,[0.25 0.5 0.75]);
    HDmeanQ_ALLenvs(pp) = mean(temp_HD_allenvs);
    JDmeanQ_ALLenvs(pp) = mean(temp_JD_allenvs);

end

% %% cdf plots:
% 
% figure
% hold on
% h(1,1) = cdfplot(HDmedQ_compenvs(:,1));
% h(1,2) = cdfplot(HDmedQ_compenvs(:,2));
% h(1,3) = cdfplot(HDmedQ_compenvs(:,3));
% set( h(1,:), 'Color', [208/247 28/247 139/247]);
% h(2,1) = cdfplot(HDmedQ_fcoopenvs(:,1));
% h(2,2) = cdfplot(HDmedQ_fcoopenvs(:,2));
% h(2,3) = cdfplot(HDmedQ_fcoopenvs(:,3));
% set( h(2,:), 'Color', [77/247 172/247 38/247]);
% % h(3,1) = cdfplot(HDmedQ_ALLenvs(:,1));
% % h(3,2) = cdfplot(HDmedQ_ALLenvs(:,2));
% % h(3,3) = cdfplot(HDmedQ_ALLenvs(:,3));
% h(3,1) = cdfplot(HDmedQ_compvsfcoopenvs(:,1));
% h(3,2) = cdfplot(HDmedQ_compvsfcoopenvs(:,2));
% h(3,3) = cdfplot(HDmedQ_compvsfcoopenvs(:,3));
% set( h(3,:), 'Color', 'k');
% set( h(:,1), 'LineStyle', '--', 'LineWidth', 3);
% set( h(:,2), 'LineStyle', '-', 'LineWidth', 3);
% set( h(:,3), 'LineStyle', '--', 'LineWidth', 3);
% title('')
% xlabel('Hamming distance between environments')
% ylabel('CDF')
% set(gca, 'FontSize', 20)
% 
% % figure
% % hold on
% % h(1,1) = cdfplot(HDmedQ_compenvs(:,1)./HDmedQ_ALLenvs(:,2));
% % h(1,2) = cdfplot(HDmedQ_compenvs(:,2)./HDmedQ_ALLenvs(:,2));
% % h(1,3) = cdfplot(HDmedQ_compenvs(:,3)./HDmedQ_ALLenvs(:,2));
% % set( h(1,:), 'Color', [208/247 28/247 139/247]);
% % h(2,1) = cdfplot(HDmedQ_fcoopenvs(:,1)./HDmedQ_ALLenvs(:,2));
% % h(2,2) = cdfplot(HDmedQ_fcoopenvs(:,2)./HDmedQ_ALLenvs(:,2));
% % h(2,3) = cdfplot(HDmedQ_fcoopenvs(:,3)./HDmedQ_ALLenvs(:,2));
% % set( h(2,:), 'Color', [77/247 172/247 38/247]);
% % h(3,1) = cdfplot(HDmedQ_ALLenvs(:,1)./HDmedQ_ALLenvs(:,2));
% % h(3,2) = cdfplot(HDmedQ_ALLenvs(:,2)./HDmedQ_ALLenvs(:,2));
% % h(3,3) = cdfplot(HDmedQ_ALLenvs(:,3)./HDmedQ_ALLenvs(:,2));
% % set( h(3,:), 'Color', 'k');
% % set( h(:,1), 'LineStyle', '--', 'LineWidth', 3);
% % set( h(:,2), 'LineStyle', '-', 'LineWidth', 3);
% % set( h(:,3), 'LineStyle', '--', 'LineWidth', 3);
% % title('')
% % xlabel('Hamming distance between environments')
% % ylabel('CDF')
% % set(gca, 'FontSize', 20)
% 
% figure
% hold on
% h(1) = cdfplot(HDmeanQ_compenvs./HDmeanQ_ALLenvs);
% h(2) = cdfplot(HDmeanQ_fcoopenvs./HDmeanQ_ALLenvs);
% % h(3) = cdfplot(HDmeanQ_ALLenvs);
% h(3) = cdfplot(HDmeanQ_compvsfcoopenvs./HDmeanQ_ALLenvs);
% set( h(1), 'Color', [208/247 28/247 139/247]);
% set( h(2), 'Color', [77/247 172/247 38/247]);
% set( h(3), 'Color', 'k');
% set( h(1), 'LineStyle', '--', 'LineWidth', 3);
% set( h(2), 'LineStyle', '--', 'LineWidth', 3);
% set( h(3), 'LineStyle', '--', 'LineWidth', 3);
% title('')
% xlabel('mean Hamming distance between environments')
% ylabel('CDF')
% set(gca, 'FontSize', 20)
% 
% 
% % figure
% % hold on
% % h(1,1) = cdfplot(JDmedQ_compenvs(:,1));
% % h(1,2) = cdfplot(JDmedQ_compenvs(:,2));
% % h(1,3) = cdfplot(JDmedQ_compenvs(:,3));
% % set( h(1,:), 'Color', [208/247 28/247 139/247]);
% % h(2,1) = cdfplot(JDmedQ_fcoopenvs(:,1));
% % h(2,2) = cdfplot(JDmedQ_fcoopenvs(:,2));
% % h(2,3) = cdfplot(JDmedQ_fcoopenvs(:,3));
% % set( h(2,:), 'Color', [77/247 172/247 38/247]);
% % h(3,1) = cdfplot(JDmedQ_ALLenvs(:,1));
% % h(3,2) = cdfplot(JDmedQ_ALLenvs(:,2));
% % h(3,3) = cdfplot(JDmedQ_ALLenvs(:,3));
% % set( h(3,:), 'Color', 'k');
% % set( h(:,1), 'LineStyle', '--', 'LineWidth', 3);
% % set( h(:,2), 'LineStyle', '-', 'LineWidth', 3);
% % set( h(:,3), 'LineStyle', '--', 'LineWidth', 3);
% % title('')
% % xlabel('Jaccard distance between environments')
% % ylabel('CDF')
% % set(gca, 'FontSize', 20)
% 
% % figure
% % hold on
% % h(1) = cdfplot(JDmeanQ_compenvs);
% % h(2) = cdfplot(JDmeanQ_fcoopenvs);
% % h(3) = cdfplot(JDmeanQ_ALLenvs);
% % set( h(1), 'Color', [208/247 28/247 139/247]);
% % set( h(1), 'Color', [77/247 172/247 38/247]);
% % set( h(1), 'Color', 'k');
% % title('')
% % xlabel('mean Jaccard distance between environments')
% % ylabel('CDF')
% % set(gca, 'FontSize', 20)
% 
% % figure
% % hold on
% % cdfplot(HDmedQ_ALLenvs(:,2));
% % cdfplot(HDmedQ_fcoopenvs(:,2));
% % cdfplot(HDmedQ_compenvs(:,2));
% % legend('all', 'fcoop', 'comp')

%%

figure
xJD = [ones(length(JDmeanQ_compenvs),1); 2*ones(length(JDmeanQ_fcoopenvs),1); 3*ones(length(JDmeanQ_compvsfcoopenvs),1)];
name_order = convertCharsToStrings([{'comp'},{'coop'},{'comp vs coop'}]);
xJD_named = categorical(xJD, 1:3, name_order);
yJD = [JDmeanQ_compenvs', JDmeanQ_fcoopenvs', JDmeanQ_compvsfcoopenvs'];
% yJD = [JDmeanQ_compenvs'./JDmeanQ_ALLenvs', JDmeanQ_fcoopenvs'./JDmeanQ_ALLenvs', JDmeanQ_compvsfcoopenvs'./JDmeanQ_ALLenvs'];
boxchart(xJD_named, yJD)
ylab = ylabel('mean Jaccard distance');
ylab.FontName = 'Times';
ylim([0,1])
set(gca,'LineWidth',3,'TickLabelInterpreter','Latex','TickLength',[.015 .015],'FontSize',20,'Box','off');
%save
filedr_savePlots = '/Users/josephinesolowiej-wedderburn/Documents/EnvProject/Plots/SI_envSimilarity/';
figstr = [filedr_savePlots, colAbr, '_meanJDwithinpairsCompCoopVs_EnvSize', num2str(nEC), '.fig'];
saveas(gcf,figstr)
% svgstr = [filedr_savePlots, colAbr, '_meanJDwithinpairsCompCoopVs_EnvSize', num2str(nEC), '.svg'];
% saveas(gcf,svgstr)
pngstr = [filedr_savePlots, colAbr, '_meanJDwithinpairsCompCoopVs_EnvSize', num2str(nEC), '.png'];
saveas(gcf,pngstr)

%check normality:
figure
normplot(JDmeanQ_compenvs)
normplot(JDmeanQ_fcoopenvs)
normplot(JDmeanQ_compvsfcoopenvs)
%
[pJD,tblJD,statsJD] = anova1(yJD, xJD);
figure
resultsJD = multcompare(statsJD);


% %% 
% 
% figure
% xHD = [ones(length(HDmeanQ_compenvs),1); 2*ones(length(HDmeanQ_fcoopenvs),1); 3*ones(length(HDmeanQ_compvsfcoopenvs),1)];
% name_order = convertCharsToStrings([{'comp'},{'coop'},{'comp vs coop'}]);
% xHD_named = categorical(xHD, 1:3, name_order);
% yHD = [HDmeanQ_compenvs', HDmeanQ_fcoopenvs', HDmeanQ_compvsfcoopenvs'];
% boxchart(xHD_named, yHD)
% ylab = ylabel('mean Hamming distance');
% ylab.FontName = 'Times';
% set(gca,'LineWidth',3,'TickLabelInterpreter','Latex','TickLength',[.015 .015],'FontSize',20,'Box','off');
% 
% %check normality:
% figure
% normplot(HDmeanQ_compenvs)
% normplot(HDmeanQ_fcoopenvs)
% normplot(HDmeanQ_compvsfcoopenvs)
% %
% [pHD,tblHD,statsHD] = anova1(yHD, xHD);
% figure
% resultsHD = multcompare(statsHD);