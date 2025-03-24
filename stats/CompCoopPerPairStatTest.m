filedr_data = '/Users/josephinesolowiej-wedderburn/Documents/EnvProject/Data_PaperFigures/';
filedr_savePlots = '/Users/josephinesolowiej-wedderburn/Documents/EnvProject/Plots/scatterTriangles/';

collection = 'AGORA';
colAbr = 'AG';
% collection = 'CarveMe';
% colAbr = 'CM';
%% %load saved data
load([filedr_data, colAbr, 'summary_data_ECs50_conc1000_100GEsx4envs.mat']);
load([filedr_data, colAbr, 'summary_data_ECs100_conc1000_100GEsx4envs.mat']);
load([filedr_data, colAbr, 'summary_data_ECs50_conc500_100GEsx4envs.mat']);
load([filedr_data, colAbr, 'summary_data_ECs100_conc500_100GEsx4envs.mat']);

%% split into bins
nbins = 10;
binsize = 10000/nbins;
whichPairs = datasample(1:10000, nbins*binsize, 'Replace', false);
binSplit = repmat(1:nbins, 1, binsize);

% %% one-way ANOVA to test numbers of comp/coop do not significantly differ across bins
% %fac cooperation:
% [pfcECs50_conc1000,tbl,stats] = anova1(summary_data_ECs50_conc1000(:,1), binSplit);
% [pfcECs100_conc1000,tbl,stats] = anova1(summary_data_ECs100_conc1000(:,1), binSplit);
% [pfcECs50_conc500,tbl,stats] = anova1(summary_data_ECs50_conc500(:,1), binSplit);
% [pfcECs100_conc500,tbl,stats] = anova1(summary_data_ECs100_conc500(:,1), binSplit);
% 
% %any cooperation:
% [pacECs50_conc1000,tbl,stats] = anova1((summary_data_ECs50_conc1000(:,1) + summary_data_ECs50_conc1000(:,10) + summary_data_ECs50_conc1000(:,13) + summary_data_ECs50_conc1000(:,16)), binSplit);
% [pacECs100_conc1000,tbl,stats] = anova1((summary_data_ECs100_conc1000(:,1) + summary_data_ECs100_conc1000(:,10) + summary_data_ECs100_conc1000(:,13) + summary_data_ECs100_conc1000(:,16)), binSplit);
% [pacECs50_conc500,tbl,stats] = anova1((summary_data_ECs50_conc500(:,1) + summary_data_ECs50_conc500(:,10) + summary_data_ECs50_conc500(:,13) + summary_data_ECs50_conc500(:,16)), binSplit);
% [pacECs100_conc500,tbl,stats] = anova1((summary_data_ECs100_conc500(:,1) + summary_data_ECs100_conc500(:,10) + summary_data_ECs100_conc500(:,13) + summary_data_ECs100_conc500(:,16)), binSplit);
% 
% %competition:
% [pcmECs50_conc1000,tbl,stats] = anova1(summary_data_ECs50_conc1000(:,4), binSplit);
% [pcmECs100_conc1000,tbl,stats] = anova1(summary_data_ECs100_conc1000(:,4), binSplit);
% [pcmECs50_conc500,tbl,stats] = anova1(summary_data_ECs50_conc500(:,4), binSplit);
% [pcmECs100_conc500,tbl,stats] = anova1(summary_data_ECs100_conc500(:,4), binSplit);

% %% one-way ANOVA to test no pairs with comp/coop does not differ significantly across bins
% AnyCoop = (summary_data_ECs100_conc1000(:,1) + summary_data_ECs100_conc1000(:,10) + summary_data_ECs100_conc1000(:,13) + summary_data_ECs100_conc1000(:,16))>0;
% AnyComp = summary_data_ECs100_conc1000(:,4)>0;
% Both = AnyCoop.*AnyComp;
% OnlyCoop = AnyCoop==1 & AnyComp==0;
% OnlyComp = AnyCoop==0 & AnyComp==1;
% neither = AnyCoop==0 & AnyComp==0;
% cat_both = ones(10000,1);
% cat_both(OnlyCoop) = 2;
% cat_both(OnlyComp) = 3;
% cat_both(neither) = 4;
% [pbothECs100_conc1000,tbl,stats] = anova1( cat_both, binSplit );

%% test the difference in proportions of comp/coop/other interactions across the different env sizes:
tot_np = 10000;
ind_summDat = 1:tot_np;
nn = 100;
pairs_nnGEs_AllCombos = ind_summDat( sum(summary_data_ECs50_conc1000(:,1:16),2)>=nn & sum(summary_data_ECs50_conc500(:,1:16),2)>=nn & sum(summary_data_ECs100_conc1000(:,1:16),2)>=nn & sum(summary_data_ECs100_conc500(:,1:16),2)>=nn);

compdata_nnGEsECs50_conc1000_vec = summary_data_ECs50_conc1000(pairs_nnGEs_AllCombos,4);
compdata_nnGEsECs100_conc1000_vec = summary_data_ECs100_conc1000(pairs_nnGEs_AllCombos,4);
anyCoopdata_nnGEsECs50_conc1000_vec = summary_data_ECs50_conc1000(pairs_nnGEs_AllCombos,1) + summary_data_ECs50_conc1000(pairs_nnGEs_AllCombos,10) + summary_data_ECs50_conc1000(pairs_nnGEs_AllCombos,13) + summary_data_ECs50_conc1000(pairs_nnGEs_AllCombos,16);
anyCoopdata_nnGEsECs100_conc1000_vec = summary_data_ECs100_conc1000(pairs_nnGEs_AllCombos,1) + summary_data_ECs100_conc1000(pairs_nnGEs_AllCombos,10) + summary_data_ECs100_conc1000(pairs_nnGEs_AllCombos,13) + summary_data_ECs100_conc1000(pairs_nnGEs_AllCombos,16);
otherdata_ECs50_conc1000_vec = (100 - (compdata_nnGEsECs50_conc1000_vec + anyCoopdata_nnGEsECs50_conc1000_vec));
otherdata_ECs100_conc1000_vec = (100 - (compdata_nnGEsECs100_conc1000_vec + anyCoopdata_nnGEsECs100_conc1000_vec));

mat1 = [compdata_nnGEsECs50_conc1000_vec, anyCoopdata_nnGEsECs50_conc1000_vec, otherdata_ECs50_conc1000_vec];
mat2 = [compdata_nnGEsECs100_conc1000_vec, anyCoopdata_nnGEsECs100_conc1000_vec, otherdata_ECs100_conc1000_vec];

% Compute observed mean difference (Euclidean distance)
observed_diff = norm(mean(mat1) - mean(mat2));

% Permutation test setup
num_permutations = 10;%10000;
perm_diffs = zeros(num_permutations, 1);

% Combine both datasets
combined = [mat1; mat2];
num_samples = size(combined, 1);

for i = 1:num_permutations
    shuffled = combined(randperm(num_samples), :); % Shuffle rows
    perm_diffs(i) = norm(mean(shuffled(1:size(mat1,1), :)) - mean(shuffled(size(mat1,1)+1:end, :)));
end

% Compute p-value
p = mean(perm_diffs >= observed_diff);
disp(['Permutation test p-value: ', num2str(p)]);

%% test the difference in proportions of fcoop/one-way/mutual interactions across the different env sizes:
tot_np = 10000;
ind_summDat = 1:tot_np;
nn = 100;
pairs_nnGEs_AllCombos = ind_summDat( sum(summary_data_ECs50_conc1000(:,1:16),2)>=nn & sum(summary_data_ECs50_conc500(:,1:16),2)>=nn & sum(summary_data_ECs100_conc1000(:,1:16),2)>=nn & sum(summary_data_ECs100_conc500(:,1:16),2)>=nn);

sumCoopdata_nnGEsECs50_conc1000_vec = summary_data_ECs50_conc1000(pairs_nnGEs_AllCombos,1) + summary_data_ECs50_conc1000(pairs_nnGEs_AllCombos,10) + summary_data_ECs50_conc1000(pairs_nnGEs_AllCombos,13) + summary_data_ECs50_conc1000(pairs_nnGEs_AllCombos,16);
sumCoopdata_nnGEsECs100_conc1000_vec = summary_data_ECs100_conc1000(pairs_nnGEs_AllCombos,1) + summary_data_ECs100_conc1000(pairs_nnGEs_AllCombos,10) + summary_data_ECs100_conc1000(pairs_nnGEs_AllCombos,13) + summary_data_ECs100_conc1000(pairs_nnGEs_AllCombos,16);

pairs_foundCoop = ind_summDat( sumCoopdata_nnGEsECs50_conc1000_vec>0 & sumCoopdata_nnGEsECs100_conc1000_vec>0 );

pfcoop_nnGEsECs50_conc1000_vec = summary_data_ECs50_conc1000(pairs_foundCoop,1)./sumCoopdata_nnGEsECs50_conc1000_vec(pairs_foundCoop);
p1way_nnGEsECs50_conc1000_vec = (summary_data_ECs50_conc1000(pairs_foundCoop,10) + summary_data_ECs50_conc1000(pairs_foundCoop,13))./sumCoopdata_nnGEsECs50_conc1000_vec(pairs_foundCoop);
pmutual_nnGEsECs50_conc1000_vec = summary_data_ECs50_conc1000(pairs_foundCoop,16)./sumCoopdata_nnGEsECs50_conc1000_vec(pairs_foundCoop);
pfcoop_nnGEsECs100_conc1000_vec = summary_data_ECs100_conc1000(pairs_foundCoop,1)./sumCoopdata_nnGEsECs100_conc1000_vec(pairs_foundCoop);
p1way_nnGEsECs100_conc1000_vec = (summary_data_ECs100_conc1000(pairs_foundCoop,10) + summary_data_ECs100_conc1000(pairs_foundCoop,13))./sumCoopdata_nnGEsECs100_conc1000_vec(pairs_foundCoop);
pmutual_nnGEsECs100_conc1000_vec = summary_data_ECs100_conc1000(pairs_foundCoop,16)./sumCoopdata_nnGEsECs100_conc1000_vec(pairs_foundCoop);

mat1 = [pfcoop_nnGEsECs50_conc1000_vec, p1way_nnGEsECs50_conc1000_vec, pmutual_nnGEsECs50_conc1000_vec];
mat2 = [pfcoop_nnGEsECs100_conc1000_vec, p1way_nnGEsECs100_conc1000_vec, pmutual_nnGEsECs100_conc1000_vec];

% Compute observed mean difference (Euclidean distance)
observed_diff = norm(mean(mat1) - mean(mat2));

% Permutation test setup
num_permutations = 10000;
perm_diffs = zeros(num_permutations, 1);

% Combine both datasets
combined = [mat1; mat2];
num_samples = size(combined, 1);

for i = 1:num_permutations
    shuffled = combined(randperm(num_samples), :); % Shuffle rows
    perm_diffs(i) = norm(mean(shuffled(1:size(mat1,1), :)) - mean(shuffled(size(mat1,1)+1:end, :)));
end

% Compute p-value
p = mean(perm_diffs >= observed_diff);
disp(['Permutation test p-value: ', num2str(p)]);