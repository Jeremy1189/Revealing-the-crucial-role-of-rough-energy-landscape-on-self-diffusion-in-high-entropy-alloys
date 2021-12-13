clc;clear;addpath(genpath(pwd)); rng('shuffle')
% close all
%% the trained Machine Learning Model
% load model_NiFeCr_all_ratio.mat model_NiFeCr_all_ratio best_model = model_NiFeCr_all_ratio;

load  best_model_EBF_NiFeCr_nn4.mat best_model_EBF_NiFeCr
best_model_EBF = best_model_EBF_NiFeCr;
% load NiFeCr_ML_FE_all_ratio.mat best_model_FE
load best_model_FE_NiFeCr_nn4.mat best_model_FE_NiFeCr
%% the standard coordinate
% read a reference structure, and the structure difference between NiFe and 
% NiFeCr is just change the types and lattice constant

%data_file_name ='NiFe55';
%per_file_path =[pwd,'\',data_file_name,'\per_NiFe'];
%per_input= read_data(per_file_path);
% [ per_input,infor_set] = read_lmp(per_file_path);
load per_input.mat per_input
vacancy_atomID=2221;% center postition
migration_atomID=vacancy_atomID+1; % which could make sure the initial Local Env is same with the trained ML  
per_index = per_input(1,:);
per_coords =  per_input(3:5,:);
vac_coord= per_coords(:,vacancy_atomID);
lattice_constant_ini=3.488;
super_size_length=10;
fix_point=10^8;
NN1=12;%fcc 1st nn
NN_ML=4;% 4NN both from migration atom and vacancy
cos_flag=1;% 1 or 2
[vac_ID_nn_set,vac_ID_count_set,~] ...
            = vac_coords_trans(vac_coord',per_coords',per_index,lattice_constant_ini,super_size_length);
cum_vac_IDs=cumsum(vac_ID_count_set);  
nn8_vac_num = cum_vac_IDs(9)-1;%5nn is the 6th number and the 1st is the vacancy 
NN=1;
NN1_atomIDs = vac_ID_nn_set(NN+1:vac_ID_count_set(NN+1)+1);
%% tranfer to new coordinate according to the Mig-vac vector

[standard_coords,standard_axis_x,standard_axis_y]= coordinates_transfer(vacancy_atomID,migration_atomID,per_coords,...
    lattice_constant_ini,super_size_length,NN_ML,cos_flag);
axis_x_set=zeros(length(NN1_atomIDs),3);
axis_y_set=zeros(length(NN1_atomIDs),3);
stand_sort_index_set=zeros(length(NN1_atomIDs),length(standard_coords));
%% update the lattice constants of NiFeCr

load optimize_lattice_constant_set0419.mat optimize_lattice_constant_set ratio
lattice_constant_set = optimize_lattice_constant_set;
%% initalization the parameters

perfect_coord=per_coords;
perfect_types=per_input(2,:);%NiFe
N=length(per_index)-1;
% iteration
iter_times=5e4;
%kmc
T=900;
Kb= 0.025852/300 ;% ev
D0=1*10^13;%HZ
L_r=size(ratio,1);
D=zeros(L_r,1);
core_number=40;
parpool('local',core_number);
parfor num_ratio=1:L_r
% for  num_ratio=1:L_r
% for  num_ratio=1:L_r
%% initialize the vacancy coordinate for each ratio

    vac_ID= vacancy_atomID;
    mig_ID=migration_atomID;
%% update the lattice constant and coordinates

    lattice_constant=lattice_constant_set(num_ratio);% lattice constant for current ratio
    Len = lattice_constant*super_size_length;% box length
    per_coords = perfect_coord ;
    per_type=perfect_types;
    per_coords =per_coords./lattice_constant_ini.*lattice_constant;
%% initialize the MSD coordinate

    R0=per_coords;
    Rt=per_coords;
%% update the types form NiFe to NiFeCr with current ratio

    num_atoms=length(per_coords);
    element_type=[1,2,3];%1 Ni, 2 Fe,3Cr
    L= length( element_type);
    num_atoms_all=zeros(L,1);
    for i= 1:L-1
        num_atoms_all(i)= round(num_atoms*ratio(num_ratio,i));
    end
    num_atoms_all(L)= num_atoms-sum(num_atoms_all); % the number of the remain element
    cumsum_set = cumsum(num_atoms_all);
    per_type(1:cumsum_set(1))=element_type(1);% Ni%
    for j = 2: L
        per_type(  (cumsum_set(j-1)+1) :cumsum_set(j))=element_type(j);
    end
    rand_index = randperm(num_atoms);
    per_type = per_type(rand_index);
%% initialization cell and array

    t=0;% inital tim
    count0=0;count60=0;count90=0;count120=0;count180=0;
    cos_value_set=zeros(iter_times,1);
    tracks= zeros(iter_times,5); %[t,vac_x,vac_y,Vac_z,MSD]
    NN1_mig_energy_set1 =zeros(iter_times,NN1);%
    NN1_FE_energy_set=zeros(iter_times,NN1);%
    mig_index_set =zeros(iter_times,1);
    K_tot_set_rng=zeros(iter_times,1);
%     per_vac_type=zeros(iter_times,1);
%     output_FE_set=zeros(iter_times,1);
    mig_type_set = zeros(iter_times,14);% [mig_type, vac_type, 12 NN1_type]
    for count=1:iter_times
%% calculating MSD

        vac_coord=per_coords(:,vac_ID);% updated the vacancy position
        Dt =abs(Rt-R0);
        index= find(Dt>0.5*Len);% more than the half of the box length
        if ~isempty(index)&&mod(count,1e4)==0
            disp(count)
        end
        Dt(index) =Len-Dt(index);
        sum_Dt=sum(Dt.^2);%
        MSD =sum(sum_Dt)./N;
        tracks(count,:)=[t;vac_coord;MSD]';
%% match the axis-y to adjust the order of input

       [vac_ID_nn_set,~,~,~,~] ...
            = update_vac_nn_kmc(vac_ID,mig_ID,per_coords',lattice_constant,super_size_length);
       %        vac_input_ID =[mig_ID;vac_ID_nn_set(1:nn5_vac_num)];% the vac_ID
%        vac_input_ID =[vac_ID;vac_ID_nn_set(1:nn8_vac_num)];% the vac_type changed with the mig_type
%        per_vac_type(count)= per_type(vac_ID);
%        input_FE = dec2bin_alloy3( per_type(vac_input_ID));
%        output_FE=  best_model_FE(input_FE');
%        output_FE_set(count)=output_FE;        
        NN1_ID_set = vac_ID_nn_set(1:NN1);        
        stand_sort_index_set=zeros(length(NN1_atomIDs),length(standard_coords));
        input_atomID_set=zeros(NN1,length(stand_sort_index_set)+1);%1 is mig
        for num_mig_NN1 =1:NN1  
          %update mig_ID
            mig_ID = NN1_ID_set(num_mig_NN1); 
            vac_mig_vector =( per_coords(:,mig_ID)-per_coords(:,vac_ID));
            norm_vac_mig_vector=vac_mig_vector./norm(vac_mig_vector);
           %match the axis x
            [new_coords,~,~]= coordinates_transfer(vac_ID,mig_ID,per_coords,...
                lattice_constant,super_size_length,NN_ML,cos_flag);
            for num_coords = 1:length(standard_coords)
                sum_diff= sum( abs(repmat(standard_coords(:,num_coords),[1,length(standard_coords)] )- new_coords),1 );
                sum_diff =round(sum_diff.*fix_point)./fix_point;
                index_p= find(sum_diff==0);
                stand_sort_index_set(num_mig_NN1,num_coords)  =index_p;% the sort order for the migration atom during the kMC
            end
              current_vac_mig_order =  stand_sort_index_set(num_mig_NN1,:);
            %update the atom_env 
            [~,~,vac_mig_sortID_nn_set,vac_mig_nn_count_set]=update_vac_nn_kmc(vac_ID,mig_ID,per_coords',lattice_constant,super_size_length);
            cum_count =cumsum(vac_mig_nn_count_set);
            NN_vac_mig_atom_ID_set=vac_mig_sortID_nn_set(1: cum_count(NN_ML));
            NN_order_ID_input= NN_vac_mig_atom_ID_set(current_vac_mig_order);
            input_atomID_set(num_mig_NN1,:)=[mig_ID ;NN_order_ID_input]';
        end
%% predicting the energy by machine learning model
%%
% 
%         Mig_IDs=input_atomID_set(:,1);
%         per=[per_index;per_type;per_coords]';
%

        input_type_set = per_type(input_atomID_set);
        input_equal_position_set= input_type_set;
        X0=dec2bin_alloy3(input_equal_position_set);
        NN1_mig_energy = best_model_EBF(X0');
        NN1_F_energy= best_model_FE_NiFeCr(X0');
%% KMC

        NN1_FE_energy_set(count,:)= NN1_F_energy;
        NN1_mig_energy_set1(count,:)=NN1_mig_energy;
        K_set= D0.*exp(-(NN1_mig_energy)./(Kb*T));
        K_tot =sum(K_set);
        K_tot_set_rng(count)=K_tot;
        cum_k_set= cumsum(K_set);
        roulette_k_set = cum_k_set./K_tot;
        r1 =rand(1);
        mig_index = find(r1-roulette_k_set <0,1);
        r2 =rand(1);
        t = t + -1/K_tot* log(r2);
%% update the Rt and vancancy coordinate

        NN1_mig_atomID = input_atomID_set(:,1);
        mig_ID =   NN1_mig_atomID(mig_index);
        mig_index_set(count)=mig_index;
        mig_coord = per_coords(:,mig_ID);
        mig_type_set(count,:)= [per_type(vac_ID),per_type(mig_ID),per_type(NN1_mig_atomID)];
        vac_coord = per_coords(:,vac_ID);
        Rt(:,mig_ID) =Rt(:,mig_ID)+(vac_coord-mig_coord);
%% jump record

        current_jump_vector=vac_coord-mig_coord;
        index_out_L= find(current_jump_vector>=0.5*Len);
        current_jump_vector(index_out_L)=current_jump_vector(index_out_L)-Len;
        index_out_R= find(current_jump_vector<-0.5*Len);
        current_jump_vector(index_out_R)=current_jump_vector(index_out_R)+Len;
        % statistic cos
        if count==1
            forward_jump_vector=current_jump_vector;
            cos_set = cosd([0,60, 90,120,180]);
        end
        cos_value_set(count) = dot(forward_jump_vector,current_jump_vector)/( norm(forward_jump_vector)*norm(current_jump_vector));
        [min_value,index_set] = min(abs(cos_value_set(count)-cos_set));
        switch index_set
            case 1
                count0=count0+1;
            case 2
                count60=count60+1;
            case 3
                count90=count90+1;
            case 4
                count120=count120+1;
            case 5
                count180=count180+1;
            otherwise
                disp("error index")
        end
        
        forward_jump_vector=current_jump_vector;
        
        % check the periodic bangdary
        Rt_coords=Rt(:, mig_ID);
        index_out= find(Rt_coords>=Len);
        if ~isempty(index_out)
            disp(['out_right',num2str(count)]);
        end
        Rt_coords(index_out)= Rt_coords(index_out)-Len;
        
        index_out= find(Rt_coords<0);
        if ~isempty(index_out)
            disp(['out_left',num2str(count)]);
        end
        Rt_coords(index_out)= Rt_coords(index_out)+Len;
        Rt(:, mig_ID)=Rt_coords;           
        %   exchange the coordinate of vacancy and migration atom
        per_coords(:,mig_ID) = vac_coord;
        per_coords(:,vac_ID)= mig_coord;
        per_type(vac_ID)=per_type(mig_ID);
%         %update vacnacy coordinate
%         vac_coord=per_coords(:,vac_ID);
        
    end
    % save the 
    t_set  =tracks(:,1);
    MSD_set  =tracks(:,5).*(N*1e-20);%1e-20 is the A^2, when transfer it to m^2
    x=t_set;
    P=polyfit(x,MSD_set,1);
    a=0:max(x)/1000:max(x);
    f = polyval(P,a);
    D(num_ratio)=P(1)./6;
    disp(['The diffusion Rate of Ni0.',num2str(round(ratio(num_ratio,1)*1e2)),...
        'Fe0.',num2str(round(ratio(num_ratio,2)*1e2)),'Cr0.',num2str(round(ratio(num_ratio,3)*1e2)),...
        ' at T=',num2str(T),' is : '])
    disp(num2str(D(num_ratio)));
    cos{num_ratio}=cos_value_set;
    tracks_set{num_ratio}=tracks;
    k_tot_avg_set{num_ratio}= K_tot_set_rng;
    count_degree(num_ratio,:)=[count0,count60,count90,count120,count180];
    NN1_mig_energy_cell{num_ratio}=NN1_mig_energy_set1;
    vac_mig_NN1_type_set{num_ratio}=mig_type_set;
    mig_index_cell{num_ratio}=mig_index_set;  
    output_FE_set_ratios{num_ratio}=NN1_FE_energy_set;
%     per_vac_type_cell{num_ratio}=per_vac_type;
end
save ML_kMC_FE_EBF_0511_nn4.mat