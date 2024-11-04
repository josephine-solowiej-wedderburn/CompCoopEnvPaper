function output = EssentialECs_lb_forAPair_NEW_BothJustGrow( Microbe1, Microbe2, nrM1, nrM2, origmodel, num_ec, num_nec, ind1, ind2, GrowTol  )
% EssentialECs_lb_forAPair_NEW_BothJustGrow( Microbe1, Microbe2, tempmodelM1ALL, tempmodelM2ALL, origmodel, num_ec, num_nec, ind1, ind2, GrowTol  )
%Note this function is a little inefficient: check the compounds that are
%'used' in net EXPORT as well as import.

%function takes in a pair of microbes and output:
%vector of essential ECs for the given pair
M1wM2nwELmat = zeros(1+num_ec,4);
M2wM1nwELmat = zeros(1+num_ec,4);

%keep a track of any catches
catches = [];   %make note if anything didn't work

%set tolerances
tolNG = 0.001;    %tolerance for no growth
tolNC = 0.01;     %tolerance for no change in growth rate

% nrM1=size(tempmodelM1ALL.lb,1); %number of reactions in modelA
% nrM2=size(tempmodelM2ALL.lb,1); %number of reactions in modelB
tnr = nrM1+nrM2;

EC_vec = 1:num_ec;

%gurobi params:
params = struct();
params.OutputFlag = 0;

%all same shared environment:
env_rhslb = -1000*1000*ones(num_ec,1); %net import is replete (microbes can uptake as much as they could desire)
%temp models already have defaulr rhsub inbuilt
% env_rhsub = Microbe1.rhs_ext_ub+Microbe2.rhs_ext_ub;  %net output is limited as 

%%M1 and M2:
%Using original parameters, set up shared environment model for M1 and M2
%and categorise the output

%%M1 st. M2 can grow above minimal threshold
    tempmodelM1wM2set = origmodel;
    tempmodelM1wM2set.A(end,nrM1+Microbe2.bmi) = 1;
    tempmodelM1wM2set.A=sparse(tempmodelM1wM2set.A);    %gurobi wants the input to be sparse
    %additional constraint that M2 grows:
    bit = 0.001;
    extrarhs = [ (GrowTol-bit) ];
    tempmodelM1wM2set.rhs(end) = extrarhs;
    %optimise for M1
    tempmodelM1wM2set.obj(Microbe1.bmi)=-1;
    outputM1stM2nw = runTempmodelMiAndMj_nWconstraints_extras( tempmodelM1wM2set, num_ec, env_rhslb, nrM1+Microbe2.bmi );
    M1wM2nwELmat(1,1) = outputM1stM2nw.RateOptMx;
    %which environmental compounds are USED in the replete environment
    M1wM2nwELmat( 1, 3 ) = sum(outputM1stM2nw.Ecompounds_used);
    M1wM2nwELmat( 2:end, 3 ) = outputM1stM2nw.Ecompounds_used;
    %which environmental compounds _might_ be used: non-zero rows in Shared env
        M1wM2whichEnvElts0 = [Microbe1.S_ext,Microbe2.S_ext] == 0;
        M1wM2rowSums = sum(M1wM2whichEnvElts0,2); %sum of each row
        M1wM2whichRowsn0 = M1wM2rowSums ~= tnr; %zero rows will have all elts == 0 (so tnr zeros)
        M1wM2UsefulECompound_vec = EC_vec(M1wM2whichRowsn0);   %vector of row numbers where compounds appear in shared environment
        M1wM2n_useful = length(M1wM2UsefulECompound_vec);
    M1wM2nwELmat(1,4) = M1wM2n_useful;
    M1wM2nwELmat((2:end),4) = M1wM2whichRowsn0;
    %which ECs are essential/limiting
    %%remove all environmental compounds that are not useful
    M1wM2new_env_rhslb = env_rhslb;
%     M1wM2whichEnvEltsARE0 = M1wM2rowSums == tnr; %zero rows will have all elts == 0 (so tnr zeros)
%     M1wM2new_env_rhslb(M1wM2whichEnvEltsARE0) = 0;
    M1wM2nw_usedindices = EC_vec(outputM1stM2nw.Ecompounds_used==1);
    %%for the used compounds, systematically remove one at a time
    for ii = 1:M1wM2nwELmat( 1, 3 )
        iout = M1wM2nw_usedindices(ii);
        %remove ii'th compound from environment input (i.e. net IMPORT = 0)
        new_env_rhslb2 = M1wM2new_env_rhslb;
        new_env_rhslb2(iout) = 0;
%         new_env_rhsub2 = env_rhsub;
        %constrain growth of M2 st. can grow above minimal threshold
        tempmodelM1wM2setNew = tempmodelM1wM2set;
        extrarhs = [ (GrowTol-bit) ];
        tempmodelM1wM2setNew.rhs(1:num_ec)=new_env_rhslb2;
        tempmodelM1wM2setNew.rhs(end)=extrarhs;
        resultM1wM2New = gurobi(tempmodelM1wM2setNew,params);
            try
                M1wM2flux1gone = abs(resultM1wM2New.objval);
                M1wM2nwELmat((iout+1),1) = M1wM2flux1gone;
%                 M1wM2nwfluxEverything = outputM1stM2nw.RateMi;
                if M1wM2flux1gone > - tolNG & M1wM2flux1gone < tolNG
                    M1wM2nwELmat((iout+1),2) = 3;
                elseif abs(M1wM2flux1gone - GrowTol) > -tolNC & abs(M1wM2flux1gone - GrowTol) < tolNC
                    M1wM2nwELmat((iout+1),2) = 0;
                elseif M1wM2flux1gone > GrowTol
                    M1wM2nwELmat((iout+1),2) = 1;
                elseif M1wM2flux1gone < GrowTol
                    M1wM2nwELmat((iout+1),2) = 2;
                else
                    disp('no category')
                end
            catch 
                M1wM2nwELmat((iout+1),2) = 3;
                cstr = "Microbe " + ind1 + " and Microbe " + ind2 + " together INFEASIBLE when trying to remove EC " + iout;
                catches = [catches, cstr];
            end
    end

%%M2 st. M1 can grow above minimal threshold
    tempmodelM2wM1set = origmodel;
    tempmodelM2wM1set.A(end,Microbe1.bmi) = 1;
    tempmodelM2wM1set.A=sparse(tempmodelM2wM1set.A);    %gurobi wants the input to be sparse
    %additional constraint that M1 grows:
    bit = 0.001;
    extrarhs = [ (GrowTol-bit) ];
    tempmodelM2wM1set.rhs(end) = extrarhs;
    %optimise for M2
    tempmodelM2wM1set.obj(nrM1+Microbe2.bmi)=-1;
    outputM2stM1nw = runTempmodelMiAndMj_nWconstraints_extras( tempmodelM2wM1set, num_ec, env_rhslb, Microbe1.bmi );
    M2wM1nwELmat(1,1) = outputM2stM1nw.RateOptMx;
    %which environmental compounds are USED in the replete environment
    M2wM1nwELmat( 1, 3 ) = sum(outputM2stM1nw.Ecompounds_used);
    M2wM1nwELmat( 2:end, 3 ) = outputM2stM1nw.Ecompounds_used;
    %which environmental compounds _might_ be used: non-zero rows in Shared env
        M2wM1whichEnvElts0 = [Microbe1.S_ext,Microbe2.S_ext] == 0;
        M2wM1rowSums = sum(M2wM1whichEnvElts0,2); %sum of each row
        M2wM1whichRowsn0 = M2wM1rowSums ~= tnr; %zero rows will have all elts == 0 (so tnr zeros)
        M2wM1UsefulECompound_vec = EC_vec(M2wM1whichRowsn0);   %vector of row numbers where compounds appear in shared environment
        M2wM1n_useful = length(M2wM1UsefulECompound_vec);
    M2wM1nwELmat(1,4) = M2wM1n_useful;
    M2wM1nwELmat((2:end),4) = M2wM1whichRowsn0;
    %which ECs are essential/limiting
    %%remove all environmental compounds that are not useful
    M2wM1new_env_rhslb = env_rhslb;
%     M2wM1whichEnvEltsARE0 = M2wM1rowSums == tnr; %zero rows will have all elts == 0 (so tnr zeros)
%     M2wM1new_env_rhslb(M2wM1whichEnvEltsARE0) = 0;
    M2wM1nw_usedindices = EC_vec(outputM2stM1nw.Ecompounds_used==1);
    %%for the used compounds, systematically remove one at a time
   for ii = 1:M2wM1nwELmat( 1, 3 )
        iout = M2wM1nw_usedindices(ii);
        %remove ii'th compound from environment input (i.e. net IMPORT = 0)
        new_env_rhslb2 = M2wM1new_env_rhslb;
        new_env_rhslb2(iout) = 0;
        %constrain growth of M1 st. above minimal threshold
        tempmodelM2wM1setNew = tempmodelM2wM1set;
        extrarhs = [ (GrowTol-bit) ];
        tempmodelM2wM1setNew.rhs(1:num_ec)=new_env_rhslb2;
        tempmodelM2wM1setNew.rhs(end)=extrarhs;
        resultM2wM1New = gurobi(tempmodelM2wM1setNew,params);
            try
                M2wM1flux1gone = abs(resultM2wM1New.objval);
                M2wM1nwELmat((iout+1),1) = M2wM1flux1gone;
                if M2wM1flux1gone > - tolNG & M2wM1flux1gone < tolNG
                    M2wM1nwELmat((iout+1),2) = 3;
                elseif abs(M2wM1flux1gone - GrowTol) > -tolNC & abs(M2wM1flux1gone - GrowTol) < tolNC
                    M2wM1nwELmat((iout+1),2) = 0;
                elseif M2wM1flux1gone > GrowTol
                    M2wM1nwELmat((iout+1),2) = 1;
                elseif M2wM1flux1gone < GrowTol
                    M2wM1nwELmat((iout+1),2) = 2;
                else
                    disp('no category')
                end
            catch
                M2wM1nwELmat((iout+1),2) = 3;
                cstr = "Microbe " + ind1 + " and Microbe " + ind2 + " together INFEASIBLE when trying to remove EC " + iout;
                catches = [catches, cstr];
            end
    end

is_essential = M1wM2nwELmat(2:end,2)==3 | M2wM1nwELmat(2:end,2)==3;

output.ess_ECs = EC_vec(is_essential);
output.catches = catches;

% %extras:
% % output.M1ELmat = M1ELmat;
% % output.M2ELmat = M2ELmat;
% output.M1wM2nwELmat = M1wM2nwELmat;
% output.M2wM1nwELmat = M2wM1nwELmat;

end