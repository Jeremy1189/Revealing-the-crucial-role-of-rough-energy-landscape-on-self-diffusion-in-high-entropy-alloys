function [vac_ID_nn_set,vac_ID_count_set,vac_mig_sortID_nn_set,vac_mig_nn_count_set,relative_vac_coords] ...
    = update_vac_nn(vac_ID,mig_ID,per_coords)
% Description:
%       This function aims to obtain the atomIDs in nn shells
% input:
%      vac_ID:size=[1,1]
%      mig_ID: size=[1,1]
%      per_coords: size=[n,3]

% output:
%       vac_mig_sortID_nn_set: 
%       vac_mig_nn_count_set
%       vac_ID_nn_set:
%       vac_ID_count_set:
% time:
%          2020/11/05 (first version)
%%
vac_coord = per_coords(vac_ID,:);
L = length(per_coords);
mig_coord = per_coords(mig_ID,:);
relative_vac_coords = per_coords - repmat(vac_coord,[L,1]);
relative_mig_coords = per_coords - repmat( mig_coord,[L,1]);
fix_point_trans=10^8;
relative_vac_coords= round(relative_vac_coords.*fix_point_trans)./fix_point_trans;
relative_mig_coords= round(relative_mig_coords.*fix_point_trans)./fix_point_trans;
% periodic boundary check
min_per_coords= min(per_coords);
% max_per_coords= max(per_coords);
lattice_constant=3.489;
supersize=8;
max_per_coords=lattice_constant*supersize*ones(1,3);
box_size =  max_per_coords-min_per_coords;
upper_boundary = round(1/2.*box_size.*fix_point_trans)./fix_point_trans;

% err_tol =10^-2;
for i= 1:length(box_size)
    vac_index_modify= find(relative_vac_coords(:,i)>upper_boundary(i));
    mig_index_modify= find( relative_mig_coords(:,i)>upper_boundary(i));
%     vac_index_modify= find(relative_vac_coords(:,i)-1/2*box_size(i)>err_tol);
%     mig_index_modify= find( relative_mig_coords(:,i)-1/2*box_size(i)>err_tol);

    relative_vac_coords(vac_index_modify,i)= relative_vac_coords(vac_index_modify,i)-box_size(i);
    relative_mig_coords(mig_index_modify,i)= relative_mig_coords(mig_index_modify,i)-box_size(i);
    % less than and equal to the box_size
%     vac_index_modify= find(relative_vac_coords(:,i)+1/2*box_size(i)<-err_tol);
%     mig_index_modify= find( relative_mig_coords(:,i)+1/2*box_size(i)<-err_tol);
    vac_index_modify= find(relative_vac_coords(:,i)<=-upper_boundary(i));
    mig_index_modify= find( relative_mig_coords(:,i)<=-upper_boundary(i));
    relative_vac_coords(vac_index_modify,i)= relative_vac_coords(vac_index_modify,i)+box_size(i);
    relative_mig_coords(mig_index_modify,i)= relative_mig_coords(mig_index_modify,i)+box_size(i);
end
% calculated the distance
vac_nn_distance = sqrt( sum(relative_vac_coords.^2,2) );
mig_nn_distance = sqrt( sum(relative_mig_coords.^2,2 ));
% % nn sort
% [sort_vac_dis,sort_vac_index] = sort(vac_nn_distance);
% [sort_mig_dis,sort_mig_index] = sort(mig_nn_distance);
% fixed
fixed_point=8;
vac_nn_distance = round(vac_nn_distance*10^fixed_point)/10^fixed_point;
mig_nn_distance = round(mig_nn_distance*10^fixed_point)/10^fixed_point;
%unique distance
%vac
unique_vac_dis_set= unique(vac_nn_distance);
L_unique = length(unique_vac_dis_set);
largest_nn=15;
vac_ID_nn_set =zeros(largest_nn,1);% output1
vac_ID_count_set= zeros(largest_nn,1);% remove vacancy itself distance 
% mig
unique_mig_dis_set= unique(mig_nn_distance);
% L_unique = length(unique_mig_dis_set); % it is euqal with the above L_unique
mig_nn_count= zeros(largest_nn,1);
mig_ID_nn_set =zeros(largest_nn,1);
count_vac=0;
count_mig=0;
unique_vac_mig_nn_ID_set  =[]; 

vac_mig_nn_count_set = zeros(largest_nn,1);
for num_unique = 1:largest_nn
    %vac
    index_vac_set=find(vac_nn_distance ==unique_vac_dis_set(num_unique+1));
    vac_ID_count_set(num_unique) = length(index_vac_set);
    vac_ID_nn_set(count_vac+1 : count_vac+ length(index_vac_set) ) = index_vac_set;%output2
    count_vac = count_vac + length(index_vac_set);
    %mig
    index_mig_set=find(mig_nn_distance  ==unique_mig_dis_set(num_unique+1));
    mig_nn_count(num_unique) = length(index_mig_set);
    mig_ID_nn_set(count_mig+1 : count_mig + length(index_mig_set) ) = index_mig_set;
    count_mig = count_mig + length(index_mig_set);
    % vac-mig    
    vac_ID_index= find(index_mig_set == vac_ID,1);
    index_mig_set(vac_ID_index) = [];
    mig_ID_index= find(index_vac_set == mig_ID,1);
    index_vac_set(mig_ID_index) = [];
    dif_set = setdiff(index_mig_set,index_vac_set);
    temp_nn_ID_set= [index_vac_set;dif_set];
    nn_ID_set= setdiff( temp_nn_ID_set, unique_vac_mig_nn_ID_set);    
    unique_vac_mig_nn_ID_set =[unique_vac_mig_nn_ID_set;nn_ID_set];% output 3
    vac_mig_nn_count_set(num_unique)= length(nn_ID_set);% %output4    
end
% output
vac_mig_sortID_nn_set=unique_vac_mig_nn_ID_set;% output 3

end

