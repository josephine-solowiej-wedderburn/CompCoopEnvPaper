function output = runTempmodelMiAndMj_nWconstraints_extras( tempmodelMiAndMj, num_ec, env_rhslb, bminW )
%function takes in the tempmodelMi from createTempmodel1Microbe_noenvrhslbs
%adds in the given env_rhslb
%and runs gurobi optimisation for the optimal growth rate

tempmodelMiAndMj.rhs(1:num_ec)=env_rhslb;

%gurobi params:
params = struct();
params.OutputFlag = 0;
%optimise:
resultTogether = gurobi(tempmodelMiAndMj,params);

lambda_optMx = abs(resultTogether.objval);
output.RateOptMx = lambda_optMx;
output.RateNoWorse = resultTogether.x(bminW);


    %which environmental compounds are present?
    %i.e. which compounds should we keep/ immediately toss out when
    %minimising
RHS = tempmodelMiAndMj.A * resultTogether.x;
EcompoundsRHS = RHS(1:num_ec);

fluxtol = 1*10^(-6);

Ecompounds_NOTused = (EcompoundsRHS<fluxtol & EcompoundsRHS>-fluxtol);
Ecompounds_used = (EcompoundsRHS>fluxtol | EcompoundsRHS<-fluxtol);
Ecompounds_usedOUT = EcompoundsRHS<-fluxtol;

    output.Ecompounds_NOTused = Ecompounds_NOTused;
    output.Ecompounds_used = Ecompounds_used;
    output.Ecompounds_usedOUT = Ecompounds_usedOUT;

    output.ECfluxMat = EcompoundsRHS;

%     %extra
%     output.rxnfluxes = resultTogether.x;
   
end