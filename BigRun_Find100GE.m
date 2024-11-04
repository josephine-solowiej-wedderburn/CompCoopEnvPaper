%Code to find 100 growth environments of each type for all the 10,000 pairs

%set-up dataset
%this is the main folder where the met model folders are found (not specific to AGORA or CarveMe):
filedr = 'sampleData/';

%choose collection and upload pairs
%% AGORA
collection = 'AGORA';     %which collection are we in
colAbr = 'AG';
fnam = [colAbr, '_tempData/'];

% %% CarveMe
% collection = 'CarveMe';     %which collection are we in
% colAbr = 'CM';
% %set up to save
% fnam = 'tempData/';

%make save folder:
eval([ 'mkdir ' fnam])

%load the pairs
%%load the pairs from BigRunFolder
eval(['load ', filedr, 'somepairs.mat']); %how to load from cluster

%search parameters
np = length(pairmatMM); %how many pairs to look at
pinds = 1:np;
%set env sizes
add_num = [ 50, 100];
num_sizes = length(add_num);
nes = 1000;   %max number of 'random' environments to search for each designated no. of additional ECs
nWantGEs = 100;
%set lb val
lb_setVals = [-1000, -500];
num_concs = length(lb_setVals);
%put these into vectors to loop over
n_param_combos = num_sizes*num_concs;
size_loop = repmat(add_num,num_concs,1);
size_loop = reshape(size_loop,1,n_param_combos);
conc_loop = repmat(lb_setVals,1,num_sizes);

%parameters for random seed
session = 0;    %number of times that this code has been run before to generate random environments
sessions = session+1;
nloops = max(pinds);

% general info about models in this collection
tM = getmetmodel_fromFile( 1, collection, filedr );
num_ec = length(tM.rhs_ext_lb);    %no. of compounds in the environment
num_nec = length(tM.rhs_int_lb);  %no. of compounds not in the environment
EC_vec = 1:num_ec;
clear tM

%tolerances
tol = 0.001;    %tolerance for categorisation of competition/cooperation
nWtol = 0.001;  %tolerance that Mi did no worse together compared to alone
g_tol = 0.1;    %this is the tolerance for growth in a replete environment for the essential compound search

% parpool;
tic
parfor ii = 1:np

%     disp(ii)

    %set up seed for each parallel loop
    s1 = RandStream.create('mrg32k3a','NumStreams',nloops*sessions,'StreamIndices',nloops*session+pinds(ii));

    %set up matrices for data
    summary_data = zeros(n_param_combos,25);
        %data to store: keep a counter of each interaction type:
        %1-16: growth categories;
        %17: no growth category
        %18: store the first growth category
        %19: store number of iterations
        %20: index of microbe 1
        %21: index of microbe 2
        %22: number of essential compounds for the pair
        %23: number of potentially useful (but not essential) compounds
        %24: env size tried
        %25: conc
        %the rows correspond to different combinations of parameters
    summary_data_f100 = zeros(n_param_combos,25);   %set up matrix to also store the first 100 data searched
    %environments:
    comp_envs = {};
    coop_envs = {};
    tog_envs = {};
    %growth rates:
    comp_GrowthRates = {};
    coop_GrowthRates = {};
    tog_GrowthRates = {};
    %any catches
    catches = [];

    %select given pair
    ind1 = pairmatMM(pinds(ii),1);
    Microbe1 = getmetmodel_fromFile( ind1, collection, filedr );
    ind2 = pairmatMM(pinds(ii),2);
    Microbe2 = getmetmodel_fromFile( ind2, collection, filedr );
    
    %calculate useful parameters:
    nrM1=size(Microbe1.S_int,2); %number of reactions in modelA
    nrM2=size(Microbe2.S_int,2); %number of reactions in modelB
    tnr = nrM1+nrM2;             %total number of reactions

    %order microbes by number of reactions
    if nrM1 < nrM2
        ind1 = pairmatMM(pinds(ii),2);
        Microbe1 = getmetmodel_fromFile( ind1, collection, filedr );
        ind2 = pairmatMM(pinds(ii),1);
        Microbe2 = getmetmodel_fromFile( ind2, collection, filedr );
        nrM1=size(Microbe1.S_int,2); %number of reactions in modelA
        nrM2=size(Microbe2.S_int,2); %number of reactions in modelB
        tnr = nrM1+nrM2;             %total number of reactions
    end
    %save indicies in summary_data (in order of 'largest' first
    summary_data(:,20) = ind1;
    summary_data(:,21) = ind2;

    %bmi optimisation indices
    bmiM1_smat = Microbe1.bmi;
    bmiM2_smat = nrM1+Microbe2.bmi;

    %upper bound on RHS set by default
    env_rhsub = Microbe1.rhs_ext_ub + Microbe2.rhs_ext_ub;

    %generate tempmodels for M1 and M2
    tempmodelM1 = createTempmetmodelMi_noenvrhslbs( Microbe1, env_rhsub, num_ec, num_nec );
    tempmodelM2 = createTempmetmodelMi_noenvrhslbs( Microbe2, env_rhsub, num_ec, num_nec );
    %now create a base shared model
    tempmodelM1M2 = createBaseTempmetmodelMiAndMj( Microbe1, Microbe2, num_ec, num_nec, env_rhsub );

    %determine default growth for the pair
    env_rhslb_def = Microbe1.rhs_ext_lb + Microbe2.rhs_ext_lb;
    outputDef = GrowPairFromTempModels_noWTog_extras( tempmodelM1, tempmodelM2, tempmodelM1M2, num_ec, env_rhslb_def, bmiM1_smat, bmiM2_smat, nWtol  );
    %save the rates
    DefGrowthRates = outputDef.RatesVec;
    %what is the category
    outputCat = CategoriseNEW_withStrongComp( outputDef.RatesVec, nWtol, ii);
    summary_data(:,18) = outputCat.cat;

    %identify essential compounds for the pair
    outputEss = EssentialECs_lb_forAPair_NEW_BothJustGrow( Microbe1, Microbe2, nrM1, nrM2, tempmodelM1M2, num_ec, num_nec, ind1, ind2, g_tol  );
    %save ess ECs
    EssentialECs = outputEss.ess_ECs ;
    essentials = EssentialECs;
    EEC_catches = outputEss.catches;
    catches = [catches,EEC_catches];
    %save number of ECs
    summary_data(:,22) = length(EssentialECs);

    %only want to add ECs that _might_ be used: non-zero rows in Shared env
        Sh_env = [Microbe1.S_ext,Microbe2.S_ext];
        whichEnvElts0 = Sh_env == 0;
        rowSums = sum(whichEnvElts0,2); %sum of each row
        %here remove the rows which are esential (by setting their row sums to tnr)
        rowSums(EssentialECs) = tnr;
        whichRowsn0 = rowSums ~= tnr; %zero rows will have all elts == 0 (so tnr zeros)
        UsefulECompound_vec = EC_vec(whichRowsn0);   %vector of row numbers where compounds appear in shared environment
        n_useful = length(UsefulECompound_vec);
    summary_data(:,23) = n_useful;
    %make sure that we aren't trying to add more ECs than sensible
    sensible_size = add_num < (n_useful - 10);
    add_num_thisPair = add_num(sensible_size);
    maxAddNumTP = max(add_num_thisPair);

    for pp = 1:n_param_combos
        numaddECs = size_loop(pp);
        lb_setVal = conc_loop(pp);
        summary_data(pp,24) = numaddECs;
        summary_data(pp,25) = lb_setVal;
        %check that we can look at envs this size
        if numaddECs <= maxAddNumTP
            %keep a count for this choice of params
            t_count = 0;
            coop_count = 0;
            comp_count = 0;
            GE_count = 0;
            temp_ntries = 0;
            %set up env
            starter_env_rhslb = zeros(num_ec,1);
            starter_env_rhslb(EssentialECs) = lb_setVal;    %always include essential compounds
            %loop through (max no.) different environments until find the
            %desired amount (or exhaust)
            for ies = 1:nes
                %change rhs lower bounds on given number of the useful environmental compounds
                bounds_to_change = datasample(s1,UsefulECompound_vec, numaddECs,'Replace',false);    %vector of the env_compound rows to change 
                NEWenv_rhslb = starter_env_rhslb;
                %for the chosen lb elements, change the bound to -1000
                NEWenv_rhslb(bounds_to_change) = lb_setVal;
                try
                    %try to grow microbes in new environment and store interaction
                    outputNEWEnv = GrowPairFromTempModels_noWTog_extras( tempmodelM1, tempmodelM2, tempmodelM1M2, num_ec, NEWenv_rhslb, bmiM1_smat, bmiM2_smat, nWtol  );
                    %what category is the new growth?
                    outputNewCat = CategoriseNEW_withStrongComp( outputNEWEnv.RatesVec, nWtol, ii);
                    new_category = outputNewCat.cat;
                    %keep track of categories:
                    summary_data(pp,new_category) = summary_data(pp,new_category) + 1;
                    temp_ntries = temp_ntries+1;  
                    %update growth counter
                    if new_category<17
                        GE_count = GE_count + 1;
                    end
                    %if at least one needs the other save
                    if new_category>=10 & new_category<=16
                            t_count = t_count+1;
                            tog_envs{pp}(t_count,:) = [ ind1, ind2, new_category, NEWenv_rhslb' ];
                            tog_GrowthRates{pp}(t_count,:) = outputNEWEnv.RatesVec;
                    end 
                    %save all cooperative environments
                    if new_category==1  %cooperation
                            coop_count = coop_count+1;
                            coop_envs{pp}(coop_count,:) = [ ind1, ind2, NEWenv_rhslb' ];
                            coop_GrowthRates{pp}(coop_count,:) = outputNEWEnv.RatesVec;
                    end 
                    %save all competitive environments
                    if new_category==4  %cooperation
                            comp_count = comp_count+1;
                            comp_envs{pp}(comp_count,:) = [ ind1, ind2, NEWenv_rhslb' ];
                            comp_GrowthRates{pp}(comp_count,:) = outputNEWEnv.RatesVec;
                    end  

                    summary_data(pp,19) = temp_ntries;
                    
                    % save counts and categories for first 100
                    if ies==100
                        summary_data_f100(pp,:) = summary_data(pp,:);
                    end

                    %if we have enough growth, end the search
                    if GE_count>=nWantGEs
                        break
                    end

                catch
 
                    summary_data(pp,17) = summary_data(pp,17)+1;  %count the categories
                    cstr = "env " + ies + " infeasible for Microbes " + ind1 + " and " + ind2;
                    catches = [catches, cstr];
                    temp_ntries = temp_ntries+1;    %don't forget to update counter of no envs tried

                    summary_data(pp,19) = temp_ntries;

                    %if make a catch on the 100th run also want to save data here
                    if ies==100
                        summary_data_f100(pp,:) = summary_data(pp,:);
                    end
                end

            end
            
            %clear counters
            t_count = [];
            coop_count = [];
            comp_count = [];
            GE_count = [];
            temp_ntries = [];
        end
    end

    %save everything:
    iistr = num2str(pinds(ii));
    ind1str = num2str(ind1);
    ind2str = num2str(ind2);
    % filename_as_string_start = [filedr, 'NewCode/BigRunFolder_Find100GE/', fnam, '/p', iistr, 'M', ind1str 'M', ind2str ];
    filename_as_string_start = [ fnam, '/p', iistr, 'M', ind1str 'M', ind2str ];
    %summary data
    filename_as_string = [filename_as_string_start, '_summary_data.mat' ];
    saveFunForParloop( filename_as_string, summary_data )
    filename_as_string = [filename_as_string_start, '_summary_data_f100.mat' ];
    saveFunForParloop( filename_as_string, summary_data_f100 )
    %growth rates
    filename_as_string = [filename_as_string_start, '_comp_GrowthRates.mat' ];
    saveFunForParloop( filename_as_string, comp_GrowthRates )
    filename_as_string = [filename_as_string_start, '_coop_GrowthRates.mat' ];
    saveFunForParloop( filename_as_string, coop_GrowthRates )
    filename_as_string = [filename_as_string_start, '_tog_GrowthRates.mat' ];
    saveFunForParloop( filename_as_string, tog_GrowthRates )
    filename_as_string = [filename_as_string_start, '_DefGrowthRates.mat' ];
    saveFunForParloop( filename_as_string, DefGrowthRates )
    %environments
    filename_as_string = [filename_as_string_start, '_comp_envs.mat' ];
    saveFunForParloop( filename_as_string, comp_envs )
    filename_as_string = [filename_as_string_start, '_coop_envs.mat' ];
    saveFunForParloop( filename_as_string, coop_envs )
    filename_as_string = [filename_as_string_start, '_tog_envs.mat' ];
    saveFunForParloop( filename_as_string, tog_envs )
    %essential compounds
    filename_as_string = [filename_as_string_start, '_essentials.mat' ];
    saveFunForParloop( filename_as_string, essentials )
    %catches
    filename_as_string = [filename_as_string_start, '_catches.mat' ];
    saveFunForParloop( filename_as_string, catches )
    
end
toc
% delete(gcp('nocreate'))

        





