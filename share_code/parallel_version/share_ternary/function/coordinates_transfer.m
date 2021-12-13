function [new_coords,standard_axis_x,standard_axis_y]= coordinates_transfer(vacancy_atomID,migration_atomID,per_coords,...
    lattice_constant,super_size_length,NN_ML,cos_flag)
%% Descrition:
% transfer the coordinates to the new coordinates with the Mig-Vac vector as x axis, the vertical vac-NN1-atom
% as the y axis, which is determine by the cos_flag, the z axis is
% determined by the cross product with the x and y unit vector
%% Input:
% vacancy_atomID:
% migration_atomID:
% per_coords: perfrct coordinates
% lattice_constant:
%super_size_length:
% NN_ML: this is the sorted atomIDs which contains the NN atoms both from
%        vacancy and migration atom, for example there are 18 NN1 atoms for Fcc vancacy and migration atom
%cos_flag: 1 or 2; there are 2 cosine values equal 0, which is vertical with theMig_Vac Vector. the 1 means the first
% finding one, and 2 means the second finding one, we may try this by
% obtain the same standard_axis_y during the kMC;
%% Output:
%new_coords:
% standard_axis_x:we could use this vector as the match set during the kMC
%standard_axis_y:we could use this vector as the match set during the kMC
%****************************************************************************************************************************
%% calculation the vacancy and migration NN information
[vac_ID_nn_set,vac_ID_count_set,vac_mig_sortID_nn_set,vac_mig_nn_count_set,relative_vac_coords] ...
    = update_vac_nn_kmc(vacancy_atomID,migration_atomID,per_coords',lattice_constant,super_size_length);
%% transfer to fractional coordinates
fractional_coords=relative_vac_coords'./lattice_constant;
NN=1;
cum_vac_NN=cumsum(vac_ID_count_set);
vac_NN1_atomID = vac_ID_nn_set(1:cum_vac_NN(NN));
NN1_fractional_coords = fractional_coords(:,vac_NN1_atomID );
% NN_ML=4;
cum_vac_mig_NN= cumsum(vac_mig_nn_count_set);
ML_trained_order_fractional_coords=fractional_coords(:,vac_mig_sortID_nn_set(1:cum_vac_mig_NN(NN_ML)));
% create a standard coordinate
standard_axis_x =fractional_coords(:,migration_atomID)-fractional_coords(:,vacancy_atomID);
standard_axis_x=standard_axis_x/norm(standard_axis_x);
standard_set = NN1_fractional_coords';% 18 is the NN1 atoms in both migration atom and vacancy
cos_sita=(standard_set*standard_axis_x)./(sqrt(sum(standard_axis_x' .^2)).* sqrt(sum(standard_set .^2,2)) );
% fix the value for equal comparing for the next 'find(cos_sita==0,1)'
fix_point=1e10;
cos_sita=round(cos_sita.*fix_point)./fix_point;
index_cos = find(cos_sita==0); % there are 2 values equal 0, but I just select the first 1, this cause the order change of the NN
standard_axis_y = standard_set(index_cos(cos_flag),:)'/norm(standard_set(index_cos(cos_flag),:));
standard_axis_z = cross(standard_axis_x,standard_axis_y);
% tranfer the coordinates
x_new = standard_axis_x'*ML_trained_order_fractional_coords./norm(standard_axis_x);% the projection length vector
y_new = standard_axis_y'*ML_trained_order_fractional_coords./norm(standard_axis_y);
z_new = standard_axis_z' *ML_trained_order_fractional_coords./norm(standard_axis_z );
new_coords =[x_new;y_new;z_new];
end
