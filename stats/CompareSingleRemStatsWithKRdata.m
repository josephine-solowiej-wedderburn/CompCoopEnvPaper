%script to see compare the number of compounds removed before the first
%switch (and what this is to) when we keep removing compounds vs. the data
%from single removal of compounds
%i.e. do switches happen faster or more to obligate than might be predicted
%if each compound has a fixed 'role' in the environment?

% %database
collection = 'AGORA';
oldCol = 'Agora';
colAbr = 'AG';
num_ec = 424;
% collection = 'CarveMe';
% oldCol = 'CarveMe';
% colAbr = 'CM';
% num_ec = 499;

npairs = 1000; %rather than looking at all possible pairs, choose a set number

%load data (from useKR100GE data)
filedrdat =  '/Users/josephinesolowiej-wedderburn/Documents/EnvProject/Data_PaperFigures/';
fnam_onlyU = [filedrdat, 'HPC2N_output/', collection, '_SingleRemoval_1to2500'];
%load pairs:
p_vec = 1:2500;
% load([filedrdat, colAbr, '10000pairs.mat']);
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
% %env tracker when removing one by one
% filename_as_string = [fnam_onlyU, '/Allps', 'EnvTrackerCoop.mat'];
% load(filename_as_string)
% filename_as_string = [fnam_onlyU, '/Allps', 'EnvTrackerComp.mat'];
% load(filename_as_string)
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

%count events:
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

%expected values if fixed:
% %startcomp:
% compEnvSize = cac_usedComp;
% p_comp2fcoop = mean(n_comp2coop./compEnvSize);
% p_comp2obl = mean((n_comp21wO+n_comp22wO)./compEnvSize);
% p_comp2oth = mean(n_comp2other_v1./compEnvSize);
% p_comp2noGrow = mean(n_comp2newEss./compEnvSize);
% p_comp2stay = mean(n_compStaycomp./compEnvSize);
% %run simulations:
% nsims = 300*50;
% compstart_firsttransition = zeros(nsims,1);
% compstart_firsttransitionCount = zeros(nsims,1);
% nmax = floor(mean(compEnvSize));
% for ii = 1:nsims
%     % if mod(ii,100)==0
%     %     disp(ii)
%     % end
%     for jj = 1:nmax
%         xx = rand;
%         if xx<=p_comp2fcoop
%             compstart_firsttransition(ii) = 1;  %1=f coop
%             compstart_firsttransitionCount(ii) = jj;
%             break
%         elseif xx>p_comp2fcoop && xx<=(p_comp2fcoop+p_comp2obl) 
%             compstart_firsttransition(ii) = 2;  %2=obligate
%             compstart_firsttransitionCount(ii) = jj;
%             break
%         elseif xx>(p_comp2fcoop+p_comp2obl) && xx<=(p_comp2fcoop+p_comp2obl+p_comp2oth) 
%             compstart_firsttransition(ii) = 5;  %5=other
%             compstart_firsttransitionCount(ii) = jj;
%             break
%         elseif xx>(p_comp2fcoop+p_comp2obl+p_comp2oth) && xx<=(p_comp2fcoop+p_comp2obl+p_comp2oth+p_comp2noGrow)
%             compstart_firsttransition(ii) = 3;  %3=no grow
%             compstart_firsttransitionCount(ii) = jj;
%             break
%         end
%     end
% end
% sum_comp2fcoop = sum(compstart_firsttransition==1);
% sum_comp2obl = sum(compstart_firsttransition==2);
% sum_comp2oth = sum(compstart_firsttransition==5);
% sum_comp2noGrow = sum(compstart_firsttransition==3);

%comp:
np = 300;
nr = 50;
compstart_firsttransition = zeros(np*nr,1);
compstart_firsttransitionCount = zeros(np*nr,1);
for ii = 1:np
    %fcoop:
    compEnvSize = cac_usedComp(ii);
    p_comp2fcoop = n_comp2coop(ii)./compEnvSize;
    p_comp2obl = (n_comp21wO(ii)+n_comp22wO(ii))./compEnvSize;
    p_comp2oth = n_comp2other_v1(ii)./compEnvSize;
    p_comp2noGrow = n_comp2newEss(ii)./compEnvSize;
    p_comp2stay = n_compStaycomp(ii)./compEnvSize;
    nmax = floor(mean(compEnvSize));
    for kk = 1:nr
        for jj = 1:nmax
            xx = rand;
            if xx<=p_comp2fcoop
                compstart_firsttransition(nr*(ii-1)+kk) = 1;  %1=f coop
                compstart_firsttransitionCount(nr*(ii-1)+kk) = jj;
                break
            elseif xx>p_comp2fcoop && xx<=(p_comp2fcoop+p_comp2obl) 
                compstart_firsttransition(nr*(ii-1)+kk) = 2;  %2=obligate
                compstart_firsttransitionCount(nr*(ii-1)+kk) = jj;
                break
            elseif xx>(p_comp2fcoop+p_comp2obl) && xx<=(p_comp2fcoop+p_comp2obl+p_comp2oth) 
                compstart_firsttransition(nr*(ii-1)+kk) = 5;  %5=other
                compstart_firsttransitionCount(nr*(ii-1)+kk) = jj;
                break
            elseif xx>(p_comp2fcoop+p_comp2obl+p_comp2oth) && xx<=(p_comp2fcoop+p_comp2obl+p_comp2oth+p_comp2noGrow)
                compstart_firsttransition(nr*(ii-1)+kk) = 3;  %3=no grow
                compstart_firsttransitionCount(nr*(ii-1)+kk) = jj;
                break
            end
        end
    end
end

% %fcoop:
% coopEnvSize = cac_usedCoop;
% p_coop2comp = mean(n_coop2comp./coopEnvSize);
% p_coop2obl = mean((n_coop21wO+n_coop22wO)./coopEnvSize);
% p_coop2oth = mean(n_coop2other_v1./coopEnvSize);
% p_coop2noGrow = mean(n_coop2newEss./coopEnvSize);
% p_coop2stay = mean(n_coopStaycoop./coopEnvSize);
% %run simulations:
% nsims = 300*50;
% coopstart_firsttransition = zeros(nsims,1);
% coopstart_firsttransitionCount = zeros(nsims,1);
% nmax = floor(mean(coopEnvSize));
% for ii = 1:nsims
%     % if mod(ii,100)==0
%     %     disp(ii)
%     % end
%     for jj = 1:nmax
%         xx = rand;
%         if xx<=p_coop2comp
%             coopstart_firsttransition(ii) = 4;  %4=comp
%             coopstart_firsttransitionCount(ii) = jj;
%             break
%         elseif xx>p_coop2comp && xx<=(p_coop2comp+p_coop2obl)
%             coopstart_firsttransition(ii) = 2;  %2=obligate
%             coopstart_firsttransitionCount(ii) = jj;
%             break
%         elseif xx>(p_coop2comp+p_coop2obl) && xx<=(p_coop2comp+p_coop2obl+p_coop2oth)
%             coopstart_firsttransition(ii) = 5;  %5=other
%             coopstart_firsttransitionCount(ii) = jj;
%             break
%         elseif xx>(p_coop2comp+p_coop2obl+p_coop2oth) && xx<=(p_coop2comp+p_coop2obl+p_coop2oth+p_coop2noGrow)
%             coopstart_firsttransition(ii) = 3;  %3=no grow
%             coopstart_firsttransitionCount(ii) = jj;
%             break
%         end
%     end
% end
% sum_coop2comp = sum(coopstart_firsttransition==4);
% sum_coop2obl = sum(coopstart_firsttransition==2);
% sum_coop2oth = sum(coopstart_firsttransition==5);
% sum_coop2noGrow = sum(coopstart_firsttransition==3);

%fcoop:
np = 300;
nr = 50;
coopstart_firsttransition = zeros(np*nr,1);
coopstart_firsttransitionCount = zeros(np*nr,1);
% nmax = floor(mean(coopEnvSize));
for ii = 1:np
    %fcoop:
    coopEnvSize = cac_usedCoop(ii);
    p_coop2comp = n_coop2comp(ii)./coopEnvSize;
    p_coop2obl = (n_coop21wO(ii)+n_coop22wO(ii))./coopEnvSize;
    p_coop2oth = n_coop2other_v1(ii)./coopEnvSize;
    p_coop2noGrow = n_coop2newEss(ii)./coopEnvSize;
    p_coop2stay = n_coopStaycoop(ii)./coopEnvSize;
    nmax = floor(mean(coopEnvSize));
    for kk = 1:nr
        for jj = 1:nmax
            xx = rand;
            if xx<=p_coop2comp
                coopstart_firsttransition(nr*(ii-1)+kk) = 4;  %4=comp
                coopstart_firsttransitionCount(nr*(ii-1)+kk) = jj;
                break
            elseif xx>p_coop2comp && xx<=(p_coop2comp+p_coop2obl)
                coopstart_firsttransition(nr*(ii-1)+kk) = 2;  %2=obligate
                coopstart_firsttransitionCount(nr*(ii-1)+kk) = jj;
                break
            elseif xx>(p_coop2comp+p_coop2obl) && xx<=(p_coop2comp+p_coop2obl+p_coop2oth)
                coopstart_firsttransition(nr*(ii-1)+kk) = 5;  %5=other
                coopstart_firsttransitionCount(nr*(ii-1)+kk) = jj;
                break
            elseif xx>(p_coop2comp+p_coop2obl+p_coop2oth) && xx<=(p_coop2comp+p_coop2obl+p_coop2oth+p_coop2noGrow)
                coopstart_firsttransition(nr*(ii-1)+kk) = 3;  %3=no grow
                coopstart_firsttransitionCount(nr*(ii-1)+kk) = jj;
                break
            end
        end
    end
end


%% load KR summary data and compare:

%from EnvProject folders:
filedrdat =  '/Users/josephinesolowiej-wedderburn/Documents/EnvProject/Data_PaperFigures/';
fnam_onlyU = [filedrdat, 'HPC2N_output/', collection, '_KRonlyUsed_1to500_more'];
%load the pairs
eval(['load ', filedrdat, colAbr, '10000pairs.mat']);%J laptop

npairs = 500;
same_num = 300;
nruns = 50;

%start environment cooperative:
interaction = 'coop';
int_ind = 1;

%which pairs have desired interaction
eval(['load ', filedrdat, collection, '_bigRun100GEsCombined_summary_data_ECs100_conc1000.mat']);%J laptop
DidThisPairInt = summary_data_ECs100_conc1000(1:500,int_ind)>0;
% nIntPs = sum(DidThisPairInt);
AllPairInd = 1:npairs;
IntPairInd = AllPairInd(DidThisPairInt);
%want to select the same number of random pairs throughout:
IntPairInd = datasample(IntPairInd, same_num, 'replace', false);
nIntPs = same_num;

%go through each pair with desired interaction and load the sequences of transitions
StartRuns = cell(1,nruns*nIntPs);
chainLengths = zeros(1,nruns*nIntPs);

%go through each pair and find the values
for ii = 1:nIntPs

    cc = IntPairInd(ii);

    ind1 = pairmatMM(cc,1);
    ind2 = pairmatMM(cc,2);

    %load environments for this pair:
    ThisPairStr = ['p' num2str(cc) 'M' num2str(ind1) 'M' num2str(ind2) ];
    try
        stringToLoad = [fnam_onlyU, '/', ThisPairStr, interaction, 'KRonlyU_Mat.mat'];
        load(stringToLoad)
        eval(['thisSet = ' interaction 'KRonlyU_Mat;'])
        % output = parloadfun( stringToLoad ); 
        % eval(['thisSet = output.' interaction 'KRonlyU_Mat;'])
%         thisSet = output.compKRonlyU_Mat;
    catch
        try
            ThisPairStr = ['p' num2str(cc) 'M' num2str(ind2) 'M' num2str(ind1) ];
            stringToLoad = [fnam_onlyU, '/', ThisPairStr, interaction, 'KRonlyU_Mat.mat'];
            load(stringToLoad)
            eval(['thisSet = ' interaction 'KRonlyU_Mat;'])
            % output = parloadfun( stringToLoad ); 
            % eval(['thisSet = output.' interaction 'KRonlyU_Mat;'])
%             thisSet = output.compKRonlyU_Mat;
        catch
            disp( ['pair ' num2str(cc) ' no comp'] );
        end
    end

    thisSet_justStates = thisSet(:,2:2:end);
    %add the starting state (in case it switches on the removal of the
    %first compound)
    thisSet_justStates = [int_ind*ones(1,nruns);thisSet_justStates];

    %note, some pairs end with '18' rather than '17': this was if the
    %removal of a compound resulted in an infaesible solution
    %for now, remove these
    thisSet_justStates(thisSet_justStates==18) = 17;

    %reduce the number of states 1-5:
    thisSet_redStates = thisSet_justStates; %1: coop; 4: comp
    thisSet_redStates(thisSet_justStates>=10 & thisSet_justStates<=16) = 2; %2: obligate
    thisSet_redStates(thisSet_justStates==17) = 3;  %3: no grow
    thisSet_redStates(thisSet_justStates==2 | thisSet_justStates==3 | thisSet_justStates==5 | thisSet_justStates==6 | thisSet_justStates==7 ) = 5;  %other: 5   

    for jj = 1:nruns
        V = thisSet_redStates(:,jj); %choose a series of run data
        short_vec = V(V>0);
        V = [short_vec];
        [~,~,c] = unique(V);
        t = diff([0;c])~= 0;
        ix = cumsum(t);
        out = [V(t),accumarray(ix(:),1)];

        %note, some pairs can survive on just the +100 'essential' compounds
        if out(end,1)~=3
            out = [out;
                    3, 1];
        end

        StartRuns(nruns*(ii-1)+jj) = {out};

        chainLengths(nruns*(ii-1)+jj) = size(out,1);

    end
end

maxChainLength = max(chainLengths);

FirstTransitionStartCoop_vec = zeros(nIntPs*nruns,1);
HowLongToFirstTransitionStartCoop_vec = zeros(nIntPs*nruns,1);

%go through each run and identify the first transition and how long it takes
for ii = 1:nIntPs*nruns

    out = cell2mat(StartRuns(ii));

    if out(1,1)~=int_ind
        FirstTransitionStartCoop_vec(ii) = out(1,1);
        HowLongToFirstTransitionStartCoop_vec(ii) = 1;
        break
    else
        FirstTransitionStartCoop_vec(ii) = out(2,1);
        HowLongToFirstTransitionStartCoop_vec(ii) = out(1,2)+1;
    end

end



%start environment competitive:
interaction = 'comp';
int_ind = 4;

%which pairs have desired interaction
eval(['load ', filedrdat, collection, '_bigRun100GEsCombined_summary_data_ECs100_conc1000.mat']);%J laptop
DidThisPairInt = summary_data_ECs100_conc1000(1:500,int_ind)>0;
% nIntPs = sum(DidThisPairInt);
AllPairInd = 1:npairs;
IntPairInd = AllPairInd(DidThisPairInt);
%want to select the same number of random pairs throughout:
IntPairInd = datasample(IntPairInd, same_num, 'replace', false);
nIntPs = same_num;

%go through each pair with desired interaction and load the sequences of transitions
StartRuns = cell(1,nruns*nIntPs);
chainLengths = zeros(1,nruns*nIntPs);

%go through each pair and find the values
for ii = 1:nIntPs

    cc = IntPairInd(ii);

    ind1 = pairmatMM(cc,1);
    ind2 = pairmatMM(cc,2);

    %load environments for this pair:
    ThisPairStr = ['p' num2str(cc) 'M' num2str(ind1) 'M' num2str(ind2) ];
    try
        stringToLoad = [fnam_onlyU, '/', ThisPairStr, interaction, 'KRonlyU_Mat.mat'];
        load(stringToLoad)
        eval(['thisSet = ' interaction 'KRonlyU_Mat;'])
        % output = parloadfun( stringToLoad ); 
        % eval(['thisSet = output.' interaction 'KRonlyU_Mat;'])
%         thisSet = output.compKRonlyU_Mat;
    catch
        try
            ThisPairStr = ['p' num2str(cc) 'M' num2str(ind2) 'M' num2str(ind1) ];
            stringToLoad = [fnam_onlyU, '/', ThisPairStr, interaction, 'KRonlyU_Mat.mat'];
            load(stringToLoad)
            eval(['thisSet = ' interaction 'KRonlyU_Mat;'])
            % output = parloadfun( stringToLoad ); 
            % eval(['thisSet = output.' interaction 'KRonlyU_Mat;'])
%             thisSet = output.compKRonlyU_Mat;
        catch
            disp( ['pair ' num2str(cc) ' no comp'] );
        end
    end

    thisSet_justStates = thisSet(:,2:2:end);
    %add the starting state (in case it switches on the removal of the
    %first compound)
    thisSet_justStates = [int_ind*ones(1,nruns);thisSet_justStates];

    %note, some pairs end with '18' rather than '17': this was if the
    %removal of a compound resulted in an infaesible solution
    %for now, remove these
    thisSet_justStates(thisSet_justStates==18) = 17;

    %reduce the number of states 1-5:
    thisSet_redStates = thisSet_justStates; %1: coop; 4: comp
    thisSet_redStates(thisSet_justStates>=10 & thisSet_justStates<=16) = 2; %2: obligate
    thisSet_redStates(thisSet_justStates==17) = 3;  %3: no grow
    thisSet_redStates(thisSet_justStates==2 | thisSet_justStates==3 | thisSet_justStates==5 | thisSet_justStates==6 | thisSet_justStates==7 ) = 5;  %other: 5   

    for jj = 1:nruns
        V = thisSet_redStates(:,jj); %choose a series of run data
        short_vec = V(V>0);
        V = [short_vec];
        [~,~,c] = unique(V);
        t = diff([0;c])~= 0;
        ix = cumsum(t);
        out = [V(t),accumarray(ix(:),1)];

        %note, some pairs can survive on just the +100 'essential' compounds
        if out(end,1)~=3
            out = [out;
                    3, 1];
        end

        StartRuns(nruns*(ii-1)+jj) = {out};

        chainLengths(nruns*(ii-1)+jj) = size(out,1);

    end
end

maxChainLength = max(chainLengths);

FirstTransitionStartComp_vec = zeros(nIntPs*nruns,1);
HowLongToFirstTransitionStartComp_vec = zeros(nIntPs*nruns,1);

%go through each run and identify the first transition and how long it takes
for ii = 1:nIntPs*nruns

    out = cell2mat(StartRuns(ii));

    if out(1,1)~=int_ind
        FirstTransitionStartComp_vec(ii) = out(1,1);
        HowLongToFirstTransitionStartComp_vec(ii) = 1;
        break
    else
        FirstTransitionStartComp_vec(ii) = out(2,1);
        HowLongToFirstTransitionStartComp_vec(ii) = out(1,2)+1;
    end

end


%% statistical tests:

filedr_savePlots = '/Users/josephinesolowiej-wedderburn/Documents/EnvProject/Plots/SI_extraPerturbation/';

FirstTransitionStartComp_vec2compare = [FirstTransitionStartComp_vec; compstart_firsttransition];
FirstTransitionStartCoop_vec2compare = [FirstTransitionStartCoop_vec; coopstart_firsttransition];
HowLongToFirstTransitionStartComp_vec2compare = [HowLongToFirstTransitionStartComp_vec; compstart_firsttransitionCount];
HowLongToFirstTransitionStartCoop_vec2compare = [HowLongToFirstTransitionStartCoop_vec; coopstart_firsttransitionCount];

% [p1,tbl,stats] = anova1(FirstTransitionStartComp_vec2compare, [ones(length(FirstTransitionStartComp_vec),1); 2*ones(length(compstart_firsttransition),1)] );
% [p2,tbl,stats] = anova1(FirstTransitionStartCoop_vec2compare, [ones(length(FirstTransitionStartCoop_vec),1); 2*ones(length(coopstart_firsttransition),1)] );
% [p3,tbl,stats] = anova1(HowLongToFirstTransitionStartComp_vec2compare, [ones(length(HowLongToFirstTransitionStartComp_vec),1); 2*ones(length(compstart_firsttransitionCount),1)] );
% [p4,tbl,stats] = anova1(HowLongToFirstTransitionStartCoop_vec2compare, [ones(length(HowLongToFirstTransitionStartCoop_vec),1); 2*ones(length(coopstart_firsttransitionCount),1)] );

% figure
% histogram(FirstTransitionStartComp_vec)
% figure
% histogram(compstart_firsttransition)
figure
data = [sum(FirstTransitionStartComp_vec==1), sum(compstart_firsttransition==1); sum(FirstTransitionStartComp_vec==2), sum(compstart_firsttransition==2);sum(FirstTransitionStartComp_vec==5), sum(compstart_firsttransition==5); sum(FirstTransitionStartComp_vec==3), sum(compstart_firsttransition==3); 0, sum(compstart_firsttransition==0)];
x = 1:5;
b = bar(x,data*100/(np*nr),'FaceColor','flat');
b(1).FaceColor = [0.25 0.25 0.25];     %grey for other
b(2).FaceColor = [0.75 0.75 0.75];     %grey for other
set(gca,'XTickLabel',{'f. coop','obligate','other','no grow','stay comp'});
ylabel('\% simulations','FontSize',24,'Interpreter','Latex'); % obviously the label should match the plot
set(gca,'LineWidth',3,'TickLabelInterpreter','Latex','TickLength',[.025 .025],'FontSize',16,'Box','off');
ylim([0,100])
l = legend('fig. 6 data', 'fixed prob');
l.FontSize = 18;
l.Interpreter = 'Latex';
l.Location = "north east";
l.Box = false;
h=gca; 
h.XAxis.TickLength = [0 0];
%save:
figstr = [filedr_savePlots, colAbr, '_BarChartFirstSwitch_f6vsfp_startComp.fig'];
saveas(gcf,figstr)
pngstr = [filedr_savePlots, colAbr, '_BarChartFirstSwitch_f6vsfp_startComp.png'];
saveas(gcf,pngstr)

% figure
% histogram(FirstTransitionStartCoop_vec)
% figure
% histogram(coopstart_firsttransition)
figure
data = [sum(FirstTransitionStartCoop_vec==4), sum(coopstart_firsttransition==4); sum(FirstTransitionStartCoop_vec==2), sum(coopstart_firsttransition==2); sum(FirstTransitionStartCoop_vec==5), sum(coopstart_firsttransition==5); sum(FirstTransitionStartCoop_vec==3), sum(coopstart_firsttransition==3); 0, sum(coopstart_firsttransition==0)];
x = 1:5;
b = bar(x,data*100/(np*nr),'FaceColor','flat');
b(1).FaceColor = [0.25 0.25 0.25];     %grey for other
b(2).FaceColor = [0.75 0.75 0.75];     %grey for other
set(gca,'XTickLabel',{'comp','obligate','other','no grow','stay coop'});
ylabel('\% simulations','FontSize',24,'Interpreter','Latex'); % obviously the label should match the plot
set(gca,'LineWidth',3,'TickLabelInterpreter','Latex','TickLength',[.025 .025],'FontSize',16,'Box','off');
ylim([0,100])
l = legend('fig. 6 data', 'fixed prob');
l.FontSize = 18;
l.Interpreter = 'Latex';
l.Location = "north east";
l.Box = false;
h=gca; 
h.XAxis.TickLength = [0 0];
%save:
figstr = [filedr_savePlots, colAbr, '_BarChartFirstSwitch_f6vsfp_startCoop.fig'];
saveas(gcf,figstr)
pngstr = [filedr_savePlots, colAbr, '_BarChartFirstSwitch_f6vsfp_startCoop.png'];
saveas(gcf,pngstr)


% figure
% histogram(HowLongToFirstTransitionStartComp_vec2compare)
% figure
% histogram(compstart_firsttransitionCount)
figure
hold on
hf6 = histogram(HowLongToFirstTransitionStartComp_vec2compare, 'Normalization', 'percentage');
hfp = histogram(compstart_firsttransitionCount, 'Normalization', 'percentage');
hf6.FaceColor = [0.25 0.25 0.25];
hfp.FaceColor = [0.75 0.75 0.75];
xlabel('\# compounds removed','FontSize',24,'Interpreter','Latex'); 
ylabel('\% simulations','FontSize',24,'Interpreter','Latex'); 
set(gca,'LineWidth',3,'TickLabelInterpreter','Latex','TickLength',[.025 .025],'FontSize',20,'Box','off');
% ylim([0,100])
l = legend('fig. 6 data', 'fixed prob');
l.FontSize = 18;
l.Interpreter = 'Latex';
l.Location = "north east";
l.Box = false;
h=gca; 
h.XAxis.TickLength = [0 0];
%save:
figstr = [filedr_savePlots, colAbr, '_Hist2FirstSwitch_f6vsfp_startComp.fig'];
saveas(gcf,figstr)
pngstr = [filedr_savePlots, colAbr, '_Hist2FirstSwitch_f6vsfp_startComp.png'];
saveas(gcf,pngstr)

% 
% figure
% histogram(HowLongToFirstTransitionStartCoop_vec2compare)
% figure
% histogram(coopstart_firsttransitionCount)
figure
hold on
hf6 = histogram(HowLongToFirstTransitionStartCoop_vec2compare, 'Normalization', 'percentage');
hfp = histogram(coopstart_firsttransitionCount, 'Normalization', 'percentage');
hf6.FaceColor = [0.25 0.25 0.25];
hfp.FaceColor = [0.75 0.75 0.75];
xlabel('\# compounds removed','FontSize',24,'Interpreter','Latex'); 
ylabel('\% simulations','FontSize',24,'Interpreter','Latex'); 
set(gca,'LineWidth',3,'TickLabelInterpreter','Latex','TickLength',[.025 .025],'FontSize',20,'Box','off');
% ylim([0,100])
l = legend('fig. 6 data', 'fixed prob');
l.FontSize = 18;
l.Interpreter = 'Latex';
l.Location = "north east";
l.Box = false;
h=gca; 
h.XAxis.TickLength = [0 0];
%save:
figstr = [filedr_savePlots, colAbr, '_Hist2FirstSwitch_f6vsfp_startCoop.fig'];
saveas(gcf,figstr)
pngstr = [filedr_savePlots, colAbr, '_Hist2FirstSwitch_f6vsfp_startCoop.png'];
saveas(gcf,pngstr)


% %statistical tests for comp vs coop
% FirstTransitionStartCompvsCoop = [FirstTransitionStartComp_vec; FirstTransitionStartCoop_vec];
% FirstTransitionStartCompvsCoop(FirstTransitionStartCompvsCoop==4) = 1;  %re-label comp as coop for fair comparison
% [p5,tbl,stats] = anova1(FirstTransitionStartCompvsCoop, [ones(length(FirstTransitionStartComp_vec),1); 2*ones(length(FirstTransitionStartCoop_vec),1)] );
% HowLongToFirstTransitionStartCompvsCoop = [HowLongToFirstTransitionStartComp_vec2compare; HowLongToFirstTransitionStartCoop_vec2compare];
% [p6,tbl,stats] = anova1(HowLongToFirstTransitionStartCompvsCoop, [ones(length(HowLongToFirstTransitionStartComp_vec2compare),1); 2*ones(length(HowLongToFirstTransitionStartCoop_vec2compare),1)] );
% 
