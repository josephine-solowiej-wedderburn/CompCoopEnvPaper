function output = runTempmodelMi_extras( tempmodelMi, num_ec, env_rhslb )
%function takes in the tempmodelMi from createTempmodel1Microbe_noenvrhslbs
%adds in the given env_rhslb
%and runs gurobi optimisation for the optimal growth rate

%also want to output the fluxes of environmental compounds (i.e. what is
%the rhs for environmental compounds?)

tempmodelMi.rhs(1:num_ec)=env_rhslb;

%gurobi params:
params = struct();
params.OutputFlag = 0;
%optimise:
resultJustMi = gurobi(tempmodelMi,params);

lambdai = abs(resultJustMi.objval);
output.RateMi = lambdai;

    %which environmental compounds are present?
    %i.e. which compounds should we keep/ immediately toss out when
    %minimising
compounds_Migrowth = tempmodelMi.A * resultJustMi.x;
Ecompounds_Migrowth = compounds_Migrowth(1:num_ec);

fluxtol = 1*10^(-6);

Ecompounds_NOTused = (Ecompounds_Migrowth<fluxtol & Ecompounds_Migrowth>-fluxtol);
Ecompounds_used = (Ecompounds_Migrowth>fluxtol | Ecompounds_Migrowth<-fluxtol);
Ecompounds_usedImport = Ecompounds_Migrowth<-fluxtol;

    output.Ecompounds_NOTused = Ecompounds_NOTused;
    output.Ecompounds_used = Ecompounds_used;
    output.Ecompounds_usedImport = Ecompounds_usedImport;

    output.ECfluxMat = Ecompounds_Migrowth;

end