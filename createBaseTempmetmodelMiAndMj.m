function tempmodelMiAndMj = createBaseTempmetmodelMiAndMj( Microbei, Microbej, num_ec, num_nec, env_rhsub )
%function takes Mi and Mj as an input and creates a tempmodel structure for
%Gurobi optimisation
%Note that as we want to manipulte the environmental rhs lbs this is just
%set to zero as a placeholder
%envrhs ubs are given as an input parameter; this is because it is the sum
%of defaults for M1 and M2 in a shared environment
%only optimise the growth of bmiOpt
%set an additional constraint on one of the microbes to be no worse than
%when grown alone; this is given by index bmiNoWorse and value
%lambda_noWorse - tol (to allow some leniency in the numerical
%optimisation)

%calculate useful parameters:
nrMi=size(Microbei.lb,1); %number of reactions in modelA
nrMj=size(Microbej.lb,1); %number of reactions in modelB
tnr = nrMi+nrMj;

%bounds on compounds (shared in environment)
env_rhslb = zeros(num_ec, 1);   %placeholder for rhs lb:
nenv_rhsb = zeros(num_nec,1);   %note: have checked that rhs_int_lb=rhs_int_ub=0
%%TO DO%% write code that gives an error message if rhs_int_lb~=rhs_int_ub~=0

%define full lb and ub for flux constraints
modelMilb = full(Microbei.lb);
modelMjlb = full(Microbej.lb);
modelMiub = full(Microbei.ub);
modelMjub = full(Microbej.ub);

%Using original parameters, set up shared environment model for M1 and M2
%and categorise the output
%Stochimetric matrix for M1 and M2;
tempmodelMiAndMj.A = zeros( (2*(num_ec + num_nec) + 1), tnr );
%environmental compounds (remember to make duplicates for the constraints
%on compound concentrations)
tempmodelMiAndMj.A(1:num_ec,1:nrMi) = Microbei.S_ext;
tempmodelMiAndMj.A(1:num_ec,(1+nrMi):(tnr)) = Microbej.S_ext;
tempmodelMiAndMj.A( (1+num_ec):(2*num_ec), : ) = tempmodelMiAndMj.A( 1:num_ec, : );
%non-environmental compounds for microbe A
tempmodelMiAndMj.A( (2*num_ec+1):(2*num_ec+num_nec), 1:nrMi ) = Microbei.S_int;
tempmodelMiAndMj.A( (2*num_ec+num_nec+1):end-1, (nrMi+1):end ) = Microbej.S_int;
%final row of stoichiometric matrix is zeros; this is space for an extra constraint that we add later such that one of the microbes grows no worse than alone
tempmodelMiAndMj.A=sparse(tempmodelMiAndMj.A);    %gurobi wants the input to be sparse
%constraints on fluxes of all reactions
tempmodelMiAndMj.lb=[modelMilb;modelMjlb];
tempmodelMiAndMj.ub=[modelMiub;modelMjub];
%constraints on concentrations (and place holder for the additional constraint on the growth of Mx)
extrarhs = 0;
tempmodelMiAndMj.rhs=full([ env_rhslb;env_rhsub;nenv_rhsb;nenv_rhsb;extrarhs ]);   %final entry is a placeholder for the no worse growth constraint
tempmodelMiAndMj.sense = [repmat('>',num_ec,1);repmat('<',num_ec,1);repmat('=',num_nec,1);repmat('=',num_nec,1);'>'];
%set objective function to optimise for M1 + M2
tempmodelMiAndMj.obj=zeros(size(tempmodelMiAndMj.A,2),1); %zero vector, length total no. of reactions
%will need to set the index of the flux to be optimised to -1 (gurobi
%minimizes)

end