function Microbe = getmetmodel_fromFile( ind, collection, filedr )
%function takes a microbe index and the collection as input and outputs the
%metabolic model
%ind should be input as a number
%collection should be input in ''
%filedr should be input in ''

% filedr = '/Users/josephinesolowiej-wedderburn/Library/CloudStorage/OneDrive-Umeåuniversitet/Documents/MATLAB/';

eval(['load ',filedr,collection,'_New/ext_int_models/model',num2str(ind),'.mat']);
Microbe = metabolic_model;
clear hemodel

end