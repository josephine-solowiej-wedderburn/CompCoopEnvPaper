%Code to keep removing ECs from environments of XX pairs (if we found the
%right environments)

%set-up dataset
%this is the main folder where the met model folders are found (not specific to AGORA or CarveMe):
filedr = ''; %HERE PUT NAME OF PROJECT FOLDER
%choose collection and upload pairs

%% what do we want to look at:
np = 2500; %no. pairs
start_env_size = 100;
conc = -1000;

%% AGORA
collection = 'AGORA';     %which collection are we in
colAbr = 'AG';

% %% CarveMe
% collection = 'CarveMe';     %which collection are we in
% colAbr = 'CM';

fnam = [filedr, 'Output/', collection, '_SingleRemoval_1to', num2str(np)];
disp('AGORA single removAL for pairs 1 to 2500, more runs')
%make save folders:
eval([ 'mkdir ' fnam])

%% load the pairs
% eval(['load ', filedr, 'NewCode/BigRunFolder/', colAbr, '10000pairs.mat']);%J laptop
% %%load the pairs from BigRunFolder
eval(['load ', filedr, colAbr, '10000pairs.mat']); %how to load from cluster

%% general set up
% general info about models in this collection
tM = getmetmodel_fromFile( 1, collection, filedr );
num_ec = length(tM.rhs_ext_lb);    %no. of compounds in the environment
num_nec = length(tM.rhs_int_lb);  %no. of compounds not in the environment
clear tM

nWtol = 0.001;

%data from BigRunFolder_Find100GE comes in format:
%set env sizes
add_num = [ 50, 100];
num_sizes = length(add_num);
% nes = 1000;   %max number of 'random' environments to search for each designated no. of additional ECs
% nWantGEs = 100;
%set lb val
lb_setVals = [-1000, -500];
num_concs = length(lb_setVals);
%put these into vectors to loop over
n_param_combos = num_sizes*num_concs;
size_loop = repmat(add_num,num_concs,1);
size_loop = reshape(size_loop,1,n_param_combos);
conc_loop = repmat(lb_setVals,1,num_sizes);
%index of desired env size and conc:
outputInds = 1:length(conc_loop);
whichOutputInd = outputInds( conc_loop==conc & size_loop==start_env_size);

%set up things to save
%keep track of whether or not they comepted/cooperated (as a check)
didTheyComp = zeros(np,1);
didTheyFCoop = zeros(np,1); 
%%when removing just one thing
%number of used compounds; maybe interesting for another run
n_used_comp = zeros(np,1);
n_used_coop = zeros(np,1);
%count the number of categories in each of the simulations
coop_newcatcount = zeros(np,17);
comp_newcatcount = zeros(np,17);
%keep a tracker of the environments as they degrade
EnvTrackerCoop = zeros(np, num_ec);
EnvTrackerComp = zeros(np, num_ec);

%% parallel runs

% parpool;
tic
parfor ii = 1:np
% for ii = 1:10
    EC_vec = 1:num_ec;

    disp(ii)

    ind1 = pairmatMM(ii,1);
    ind2 = pairmatMM(ii,2);

    %select given pair
    Microbe1 = getmetmodel_fromFile( ind1, collection, filedr );
    Microbe2 = getmetmodel_fromFile( ind2, collection, filedr );

    %load environments for this pair:
    ThisPairStr = ['p' num2str(ii) 'M' num2str(ind1) 'M' num2str(ind2) ];
    try
        % stringToLoad = [filedr, 'NewCode/BigRunFolder_Find100GE/HPC2N_output/', collection, 'Try1_2500seed1_fromScratch/', ThisPairStr, '_coop_envs.mat'];
        stringToLoad = [filedr, 'Output/', collection, 'Try1_2500seed1_fromScratch/', ThisPairStr, '_coop_envs.mat'];
        coop_envs = parloadfun( stringToLoad );        
        % stringToLoad = [filedr, 'NewCode/BigRunFolder_Find100GE/HPC2N_output/', collection, 'Try1_2500seed1_fromScratch/', ThisPairStr, '_comp_envs.mat'];
        stringToLoad = [filedr, 'Output/', collection, 'Try1_2500seed1_fromScratch/', ThisPairStr, '_comp_envs.mat'];
        comp_envs = parloadfun( stringToLoad );
        % stringToLoad = [filedr, 'NewCode/BigRunFolder_Find100GE/HPC2N_output/', collection, 'Try1_2500seed1_fromScratch/', ThisPairStr, '_essentials.mat'];
        stringToLoad = [filedr, 'Output/', collection, 'Try1_2500seed1_fromScratch/', ThisPairStr, '_essentials.mat'];
        output = parloadfun( stringToLoad );
        essentials = output.essentials;
        catch
        try
            ThisPairStr = ['p' num2str(ii) 'M' num2str(ind2) 'M' num2str(ind1) ];
            % stringToLoad = [filedr, 'NewCode/BigRunFolder_Find100GE/HPC2N_output/', collection, 'Try1_2500seed1_fromScratch/', ThisPairStr, '_coop_envs.mat'];
            stringToLoad = [filedr, 'Output/', collection, 'Try1_2500seed1_fromScratch/', ThisPairStr, '_coop_envs.mat'];
            coop_envs = parloadfun( stringToLoad );        
            % stringToLoad = [filedr, 'NewCode/BigRunFolder_Find100GE/HPC2N_output/', collection, 'Try1_2500seed1_fromScratch/', ThisPairStr, '_comp_envs.mat'];
            stringToLoad = [filedr, 'Output/', collection, 'Try1_2500seed1_fromScratch/', ThisPairStr, '_comp_envs.mat'];
            comp_envs = parloadfun( stringToLoad );
            % stringToLoad = [filedr, 'NewCode/BigRunFolder_Find100GE/HPC2N_output/', collection, 'Try1_2500seed1_fromScratch/', ThisPairStr, '_essentials.mat'];
            stringToLoad = [filedr, 'Output/', collection, 'Try1_2500seed1_fromScratch/', ThisPairStr, '_essentials.mat'];
            output = parloadfun( stringToLoad );
            essentials = output.essentials;
        catch
            disp( ['pair ' num2str(ii) ' missing'] );
        end
    end

    comp_env = comp_envs.comp_envs;
    coop_env = coop_envs.coop_envs;

    %shared rhs ub:
    env_rhsub = Microbe1.rhs_ext_ub + Microbe2.rhs_ext_ub;

    %calculate useful parameters:
    nrM1=size(Microbe1.S_int,2); %number of reactions in modelA
    nrM2=size(Microbe2.S_int,2); %number of reactions in modelB
    tnr = nrM1+nrM2;             %total number of reactions

    %bmi optimisation indices
    bmiM1_smat = Microbe1.bmi;
    bmiM2_smat = nrM1+Microbe2.bmi;

    %set up temp models
    %generate tempmodels for M1 and M2
    tempmodelM1 = createTempmetmodelMi_noenvrhslbs( Microbe1, env_rhsub, num_ec, num_nec );
    tempmodelM2 = createTempmetmodelMi_noenvrhslbs( Microbe2, env_rhsub, num_ec, num_nec );
    %now create a base shared model
    tempmodelM1M2 = createBaseTempmetmodelMiAndMj( Microbe1, Microbe2, num_ec, num_nec, env_rhsub );

%%%COOPERATION%%%
    %if we have a coop environment, find it and perturb from here
    if length(coop_env)>=whichOutputInd
    if isempty(cell2mat(coop_env(whichOutputInd)))==0
        
        didTheyFCoop(ii) = 1;
        
        coop_env = cell2mat(coop_env(whichOutputInd));
        startCoopEnv = coop_env(1,3:end)';           

        %determine which compounds are used
        outputStartCoop = GrowPairFromTempModels_noWTog_extras( tempmodelM1, tempmodelM2, tempmodelM1M2, num_ec, startCoopEnv, bmiM1_smat, bmiM2_smat, nWtol  );
        CoopUsedECstart = outputStartCoop.UsedECs;

        %remove each used compound individually (with replacement) in this environment to see what would happen
        %identify the used compounds
        n_used_coop(ii) = sum(CoopUsedECstart);
        %essential compounds are always essential
%         coop_used_index = EC_vec(CoopUsedECstart);
        CoopUsedEC_startNotEss = CoopUsedECstart;
        CoopUsedEC_startNotEss(essentials) = 0;
        coop_used_index_NotEss = EC_vec(CoopUsedEC_startNotEss);

        %run through each used not essential compound in coop env; 
        %note this could have more essential compounds than in the replete environment as there are now fewer substitutions possible
        Coop_tempNCC = zeros(1,17);
        Temp_ETC = zeros(1, num_ec);
        ncoop = length(coop_used_index_NotEss);
            for cc = 1:ncoop
                c_remove = coop_used_index_NotEss(cc);
                try
                new_env_rhslb = startCoopEnv;
                new_env_rhslb(c_remove) = 0;
                output_coopNew = GrowPairFromTempModels_noWTog_extras( tempmodelM1, tempmodelM2, tempmodelM1M2, num_ec, new_env_rhslb, bmiM1_smat, bmiM2_smat, nWtol  );
                outputCat = CategoriseNEW_withStrongComp( output_coopNew.RatesVec, nWtol, cc);
                newCat = outputCat.cat;
                Coop_tempNCC(newCat) = Coop_tempNCC(newCat) + 1;
                Temp_ETC(c_remove) = newCat;
                catch
                    Coop_tempNCC(17) = Coop_tempNCC(17) + 1;
                    Temp_ETC(c_remove) = 17;
                    cstr = "env " + c_remove + " infeasible for Microbes " + ind1 + " and " + ind2;
                    disp(cstr)
                end
            end
            coop_newcatcount(ii, :) = Coop_tempNCC;
            EnvTrackerCoop(ii, :) = Temp_ETC;

    end
    end

    %%%COMPETITION%%%
    %if we have a comp environment, find it and perturb from here
    if length(comp_env)>=whichOutputInd %isempty(comp_env)==0
    if isempty(cell2mat(comp_env(whichOutputInd)))==0
        
        didTheyComp(ii) = 1;
        
        comp_env = cell2mat(comp_env(whichOutputInd));
        startCompEnv = comp_env(1,3:end)';

        %determine which compounds are used
        outputStartComp = GrowPairFromTempModels_noWTog_extras( tempmodelM1, tempmodelM2, tempmodelM1M2, num_ec, startCompEnv, bmiM1_smat, bmiM2_smat, nWtol  );
        CompUsedECstart = outputStartComp.UsedECs;

        %remove each used compound individually (with replacement) in this environment to see what would happen
        %identify the used compounds
        n_used_comp(ii) = sum(CompUsedECstart);
        %essential compounds are always essential
        CompUsedEC_startNotEss = CompUsedECstart;
        CompUsedEC_startNotEss(essentials) = 0;
        comp_used_index_NotEss = EC_vec(CompUsedEC_startNotEss);

        %run through each used not essential compound in coop env; 
        %note this could have more essential compounds than in the replete environment as there are now fewer substitutions possible
        Comp_tempNCC = zeros(1,17);
        Temp_ETC = zeros(1, num_ec);
        ncomp = length(comp_used_index_NotEss);
            for cc = 1:ncomp
                c_remove = comp_used_index_NotEss(cc);
                try
                new_env_rhslb = startCompEnv;
                new_env_rhslb(c_remove) = 0;
                output_compNew = GrowPairFromTempModels_noWTog_extras( tempmodelM1, tempmodelM2, tempmodelM1M2, num_ec, new_env_rhslb, bmiM1_smat, bmiM2_smat, nWtol  );
                outputCat = CategoriseNEW_withStrongComp( output_compNew.RatesVec, nWtol, cc);
                newCat = outputCat.cat;
                Comp_tempNCC(newCat) = Comp_tempNCC(newCat) + 1;
                Temp_ETC(c_remove) = newCat;
                catch
                    Comp_tempNCC(17) = Comp_tempNCC(17) + 1;
                    Temp_ETC(c_remove) = 17;
                    cstr = "env " + c_remove + " infeasible for Microbes " + ind1 + " and " + ind2;
                    disp(cstr)
                end
            end
            comp_newcatcount(ii, :) = Comp_tempNCC;
            EnvTrackerComp(ii, :) = Temp_ETC;

    end
    end

end
toc

% delete(gcp('nocreate'))

%% save output

%save number compounds used
filename_as_string = [fnam, '/Allps', 'n_used_coop.mat'];
saveFunForParloop( filename_as_string, n_used_coop )            
filename_as_string = [fnam, '/Allps', 'n_used_comp.mat'];
saveFunForParloop( filename_as_string, n_used_comp )  
%save categories when removing one by one
filename_as_string = [fnam, '/Allps', 'coop_newcatcount.mat'];
saveFunForParloop( filename_as_string, coop_newcatcount )            
filename_as_string = [fnam, '/Allps', 'comp_newcatcount.mat'];
saveFunForParloop( filename_as_string, comp_newcatcount )  
%save env tracker when removing one by one
filename_as_string = [fnam, '/Allps', 'EnvTrackerCoop.mat'];
saveFunForParloop( filename_as_string, EnvTrackerCoop )            
filename_as_string = [fnam, '/Allps', 'EnvTrackerComp.mat'];
saveFunForParloop( filename_as_string, EnvTrackerComp ) 
%save did they comp/coop
filename_as_string = [fnam, '/Allps', 'didTheyFCoop.mat'];
saveFunForParloop( filename_as_string, didTheyFCoop )            
filename_as_string = [fnam, '/Allps', 'didTheyComp.mat'];
saveFunForParloop( filename_as_string, didTheyComp )  
