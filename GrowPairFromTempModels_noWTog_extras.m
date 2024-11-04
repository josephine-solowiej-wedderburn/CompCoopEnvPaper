function output = GrowPairFromTempModels_noWTog_extras( tempmodelM1, tempmodelM2, tempmodelM1M2, num_ec, env_rhslb, bmiM1_smat, bmiM2_smat, nWtol  )

    output_JustM1 = runTempmodelMi_extras( tempmodelM1, num_ec, env_rhslb );
    output_JustM2 = runTempmodelMi_extras( tempmodelM2, num_ec, env_rhslb );  
    % optimise M1 such that M2 does no worse:
    tempmodelM1wM2nW = tempmodelM1M2;
    lambda_M2nW = output_JustM2.RateMi;
    tempmodelM1wM2nW.A(end,bmiM2_smat) = 1;
    tempmodelM1wM2nW.A = sparse(tempmodelM1wM2nW.A);
    tempmodelM1wM2nW.rhs(end) = lambda_M2nW-nWtol;
    tempmodelM1wM2nW.obj(bmiM1_smat)=-1;
    output_M1wM2nW = runTempmodelMiAndMj_nWconstraints_extras( tempmodelM1wM2nW, num_ec, env_rhslb, bmiM2_smat );  
    % optimise M1 such that M2 does no worse:
    tempmodelM2wM1nW = tempmodelM1M2;
    lambda_M1nW = output_JustM1.RateMi;
    tempmodelM2wM1nW.A(end,bmiM1_smat) = 1;
    tempmodelM2wM1nW.A = sparse(tempmodelM2wM1nW.A);
    tempmodelM2wM1nW.rhs(end) = lambda_M1nW-nWtol;
    tempmodelM2wM1nW.obj(bmiM2_smat) = -1;
    output_M2wM1nW = runTempmodelMiAndMj_nWconstraints_extras( tempmodelM2wM1nW, num_ec, env_rhslb, bmiM1_smat );
    
    output.RatesVec = [ output_JustM1.RateMi, output_M1wM2nW.RateOptMx, output_M2wM1nW.RateNoWorse, output_JustM2.RateMi, output_M1wM2nW.RateNoWorse,output_M2wM1nW.RateOptMx ];   %row vector of the growth rates 
    
    %used here means that it is net IMPORTed
    UsedECs = output_JustM1.Ecompounds_usedImport + output_JustM2.Ecompounds_usedImport + output_M1wM2nW.Ecompounds_usedOUT + output_M2wM1nW.Ecompounds_usedOUT;
    UsedECs = UsedECs>0;

    output.UsedECs = UsedECs;

    %maybe want to add more?
    %which ECs were not used, and fluxes?
    rhs_Efluxes_M1wM2nW = output_M1wM2nW.ECfluxMat;
    rhs_Efluxes_M2wM1nW = output_M2wM1nW.ECfluxMat;
    meanPairRHSenvFlux = (rhs_Efluxes_M1wM2nW + rhs_Efluxes_M2wM1nW)/2;

    output.meanPairRHSenvFlux = meanPairRHSenvFlux;


%     %more extras
%     UsedECsINorOUT = output_JustM1.Ecompounds_used + output_JustM2.Ecompounds_used + output_M1wM2nW.Ecompounds_used + output_M2wM1nW.Ecompounds_used;
%     UsedECsINorOUT = UsedECsINorOUT>0;
%     output.UsedECsINorOUT = UsedECsINorOUT;
% 
%     output.rxnfluxes_mat = [output_M1wM2nW.rxnfluxes, output_M2wM1nW.rxnfluxes];

end