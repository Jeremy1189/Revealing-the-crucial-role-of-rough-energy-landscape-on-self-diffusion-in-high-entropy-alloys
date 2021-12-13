function relative_coordinate = periodic_boundary_condition(central_coords,perfect_supercell_coords,lattice_constant,supercell_dublicate_size)
 L_perfect_supercell_coords = length(perfect_supercell_coords);
relative_coordinate = repmat(central_coords,  L_perfect_supercell_coords ,1) -  perfect_supercell_coords;
% periodic check and fix
[L_row,L_col]= size(relative_coordinate);
count=0;
for num_col=1:L_col
    for num_row = 1:L_row
        if relative_coordinate(num_row,num_col)>lattice_constant*supercell_dublicate_size/2
            relative_coordinate(num_row,num_col)= relative_coordinate(num_row,num_col)-lattice_constant*supercell_dublicate_size;
            count=count+1;
%         elseif relative_coordinate(num_row,num_col)<=-lattice_constant*supercell_dublicate_size/2
        elseif relative_coordinate(num_row,num_col)<-lattice_constant*supercell_dublicate_size/2
            relative_coordinate(num_row,num_col)= relative_coordinate(num_row,num_col)+lattice_constant*supercell_dublicate_size;
            count=count+1;
        elseif abs(relative_coordinate(num_row,num_col))==lattice_constant*supercell_dublicate_size/2
             if relative_coordinate(num_row,num_col)<0
%                  disp('relative_coordinate(num_row,num_col)=-lattice_constant*supercell_dublicate_size/2') 
                 
             else
%                  disp('relative_coordinate(num_row,num_col)=lattice_constant*supercell_dublicate_size/2') 
                  
             end
        else
        end
    end    
end
disp(count)

end