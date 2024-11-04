function tempmodelMi = createTempmetmodelMi_noenvrhslbs( Microbei, env_rhsub, num_ec, num_nec )
%function takes Microbe1 as an input and creates a tempmodel structure for
%Gurobi optimisation; the met model is input in the new structure with
%separated environmental and internal compartments
%Note that as we want to manipulte the environmental rhs lbs this is just
%set to zero as a placeholder
%envrhs ubs are given as an input parameter; this is because it is the sum
%of defaults for M1 and M2 in a shared environment

%placeholder for rhs lb:
env_rhslb = zeros(num_ec, 1);

% %note the Microbe comes with a rhs internal lb and ub; these should both be
% %zero
% %check so that we can set to zero:
% if sum(Microbei.rhs_int_lb~=0)>0
%     disp('non-zero rhs internal lower bound')
% end
% if sum(Microbei.rhs_int_ub~=0)>0
%     disp('non-zero rhs internal upper bound')
% end

%create model
    %Stochimetric matrix - this only has the reactions of M1;
    %still want to separate environmental and not compounds (used for bounds
    %on concentrations)
    tempmodelMi.A = [Microbei.S_ext;Microbei.S_ext;Microbei.S_int];
    tempmodelMi.A=sparse(tempmodelMi.A);    %gurobi wants the input to be sparse
    %constraints on fluxes of M1 reactions
    tempmodelMi.lb=full(Microbei.lb);
    tempmodelMi.ub=full(Microbei.ub);
    %constraints on concentrations - this is where we input the combined environment
    tempmodelMi.rhs=full([ env_rhslb;env_rhsub;Microbei.rhs_int_lb]);
    tempmodelMi.sense = [repmat('>',num_ec,1);repmat('<',num_ec,1);repmat('=',num_nec,1)];
    %set objective function for M1
    tempmodelMi.obj=zeros(size(tempmodelMi.A,2),1); %zero vector, length total no. of reactions
    tempmodelMi.obj(Microbei.bmi)=-1;  %optimising the flux through bmi reaction of M1

end