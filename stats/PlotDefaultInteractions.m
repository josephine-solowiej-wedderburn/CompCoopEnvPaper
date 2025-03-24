filedr_data = '/Users/josephinesolowiej-wedderburn/Documents/EnvProject/Data_PaperFigures/';

%% %load AG data
collection = 'AGORA';
colAbr = 'AG';
% %load bigRun100GEs data 
% eval(['load ' collection '_bigRun100GEsCombined_summary_data_ECs100_conc1000.mat']);
%load data for 10000 pairs with 100 GEs in all 4 types
eval(['load ' filedr_data, colAbr 'summary_data_ECs100_conc1000_100GEsx4envs.mat']);
summary_dataAG = summary_data_ECs100_conc1000;
%% %load CM data
collection = 'CarveMe';
colAbr = 'CM';
% %load bigRun100GEs data 
% eval(['load ' collection '_bigRun100GEsCombined_summary_data_ECs100_conc1000.mat']);
%load data for 10000 pairs with 100 GEs in all 4 types
eval(['load ' filedr_data, colAbr 'summary_data_ECs100_conc1000_100GEsx4envs.mat']);
summary_dataCM = summary_data_ECs100_conc1000;
%% what are the default interactions
n_fcoopAG = sum(summary_dataAG(:,18)==1);
n_compAG = sum(summary_dataAG(:,18)==4);
n_neutAG = sum(summary_dataAG(:,18)==7);
n_commensalAG = sum(summary_dataAG(:,18)==2 | summary_dataAG(:,18)==3);
if (n_fcoopAG + n_compAG + n_neutAG + n_commensalAG)~=10000
    disp('missing something')
end

n_fcoopCM = sum(summary_dataCM(:,18)==1);
n_compCM = sum(summary_dataCM(:,18)==4);
n_neutCM = sum(summary_dataCM(:,18)==7);
n_commensalCM = sum(summary_dataCM(:,18)==2 | summary_dataCM(:,18)==3);
if (n_fcoopCM + n_compCM + n_neutCM + n_commensalCM)~=10000
    disp('missing something')
end

%% split into bins and determine the variance

nbins = 10;
binsize = 10000/nbins;
whichPairs = datasample(1:10000, nbins*binsize, 'Replace', false);
binSplit = repmat(1:nbins, 1, binsize);

n_fcoopAG_pb = zeros(nbins,1);
n_compAG_pb = zeros(nbins,1);
n_neutAG_pb = zeros(nbins,1);
n_commensalAG_pb = zeros(nbins,1);
n_fcoopCM_pb = zeros(nbins,1);
n_compCM_pb = zeros(nbins,1);
n_neutCM_pb = zeros(nbins,1);
n_commensalCM_pb = zeros(nbins,1);
for ii = 1:nbins
    ptb = whichPairs(binSplit==ii);
    %AG
    n_fcoopAG_pb(ii) = sum(summary_dataAG(ptb,18)==1);
    n_compAG_pb(ii) = sum(summary_dataAG(ptb,18)==4);
    n_neutAG_pb(ii) = sum(summary_dataAG(ptb,18)==7);
    n_commensalAG_pb(ii) = sum(summary_dataAG(ptb,18)==2 | summary_dataAG(ptb,18)==3);
    %CM
    n_fcoopCM_pb(ii) = sum(summary_dataCM(ptb,18)==1);
    n_compCM_pb(ii) = sum(summary_dataCM(ptb,18)==4);
    n_neutCM_pb(ii) = sum(summary_dataCM(ptb,18)==7);
    n_commensalCM_pb(ii) = sum(summary_dataCM(ptb,18)==2 | summary_dataCM(ptb,18)==3);
end

%mean
Mean_fcoopAG_pb = mean(n_fcoopAG_pb./binsize);
Mean_compAG_pb = mean(n_compAG_pb./binsize);
Mean_neutAG_pb = mean(n_neutAG_pb./binsize);
Mean_commensalAG_pb = mean(n_commensalAG_pb./binsize);
Mean_fcoopCM_pb = mean(n_fcoopCM_pb./binsize);
Mean_compCM_pb = mean(n_compCM_pb./binsize);
Mean_neutCM_pb = mean(n_neutCM_pb./binsize);
Mean_commensalCM_pb = mean(n_commensalCM_pb./binsize);
%commensal and neutral are 'other'
Mean_othAG_pb = mean((n_neutAG_pb+n_commensalAG_pb)./binsize);
Mean_othCM_pb = mean((n_neutCM_pb+n_commensalCM_pb)./binsize);
%variance
Std_fcoopAG_pb = std(n_fcoopAG_pb./binsize);
Std_compAG_pb = std(n_compAG_pb./binsize);
Std_neutAG_pb = std(n_neutAG_pb./binsize);
Std_commensalAG_pb = std(n_commensalAG_pb./binsize);
Std_fcoopCM_pb = std(n_fcoopCM_pb./binsize);
Std_compCM_pb = std(n_compCM_pb./binsize);
Std_neutCM_pb = std(n_neutCM_pb./binsize);
Std_commensalCM_pb = std(n_commensalCM_pb./binsize);
%commensal and neutral are 'other'
Std_othAG_pb = std((n_neutAG_pb+n_commensalAG_pb)./binsize);
Std_othCM_pb = std((n_neutCM_pb+n_commensalCM_pb)./binsize);

% %one-way ANOVA to test proportions do not significantly differ across
% bins
[p,tbl,stats] = anova1(summary_dataAG(:,18), binSplit);
[p,tbl,stats] = anova1(summary_dataCM(:,18), binSplit);

%% plot



% figure
% x = [1,2];
% data = [ Mean_compAG_pb, Mean_compCM_pb; Mean_fcoopAG_pb, Mean_fcoopCM_pb; Mean_commensalAG_pb, Mean_commensalCM_pb; Mean_neutAG_pb, Mean_neutCM_pb ]';
% errhigh = [ Std_compAG_pb, Std_compCM_pb; Std_fcoopAG_pb, Std_fcoopCM_pb; Std_commensalAG_pb, Std_commensalCM_pb; Std_neutAG_pb, Std_neutCM_pb ]';
% b = bar(x,data*100,'FaceColor','flat'); 
% b(1).FaceColor = [208/255 28/255 139/255];      %comp (hot pink)
% b(2).FaceColor = [77/255 172/255 38/255];       %fac coop (green)
% b(3).FaceColor = [184/255 225/255 134/255];     %commensal (light green)
% b(4).FaceColor = [146/255 197/255 222/255];     %=/= (blue)
% 
% set(gca,'XTickLabel',{'AGORA','CarveMe'});
% hold on
% y=data*100;
% ngroups = size(y, 1);
% nbars = size(y, 2);
% % Calculating the width for each bar group
% groupwidth = min(0.8, nbars/(nbars + 1.5));
% for i = 1:nbars
%     x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%     er = errorbar(x, y(:,i), errhigh(:,i)*100, '.');
%     er.LineWidth = 3;
%     er.Color = [0 0 0];
% end
% 
% % xlabel('$p_i$','FontSize',24,'Interpreter','Latex'); % obviously the label should match the plot
% ylabel('percentage of pairs','FontSize',24,'Interpreter','Latex'); % obviously the label should match the plot
% set(gca,'LineWidth',3,'TickLabelInterpreter','Latex','TickLength',[.025 .025],'FontSize',18,'Box','off');
% ylim([0,100])
% 
% l = legend('Competition', 'Cooperation', 'Commensalism', 'Neutral');
% l.FontSize = 18;
% l.Interpreter = 'Latex';
% l.Location = "north";
% l.Box = false;


%comp/ coop / other
figure
x = [1,2];
data = [ Mean_compAG_pb, Mean_compCM_pb; Mean_fcoopAG_pb, Mean_fcoopCM_pb; Mean_othAG_pb, Mean_othCM_pb ]';
errhigh = [ Std_compAG_pb, Std_compCM_pb; Std_fcoopAG_pb, Std_fcoopCM_pb; Std_othAG_pb, Std_othCM_pb ]';
b = bar(x,data*100,'FaceColor','flat'); 
b(1).FaceColor = [208/255 28/255 139/255];      %comp (hot pink)
% b(2).FaceColor = [77/255 172/255 38/255];       %fac coop (green)
b(2).FaceColor = [184/247 225/247 134/247];     %fac coop (light green)
b(3).FaceColor = [0.75 0.75 0.75];     %grey for other

set(gca,'XTickLabel',{'AGORA','CarveMe'});
hold on
y=data*100;
ngroups = size(y, 1);
nbars = size(y, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(x, y(:,i), errhigh(:,i)*100, '.');
    er.LineWidth = 3;
    er.Color = [0 0 0];
end

% xlabel('$p_i$','FontSize',24,'Interpreter','Latex'); % obviously the label should match the plot
ylabel('percentage of pairs','FontSize',24,'Interpreter','Latex'); % obviously the label should match the plot
set(gca,'LineWidth',3,'TickLabelInterpreter','Latex','TickLength',[.025 .025],'FontSize',20,'Box','off');
ylim([0,100])

l = legend('Competition', 'Cooperation', 'Other');
l.FontSize = 18;
l.Interpreter = 'Latex';
l.Location = "north east";
l.Box = false;

h=gca; 
h.XAxis.TickLength = [0 0];