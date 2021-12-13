function read_data = read_data(file_path)
%current_path = C:\Users\86135\Desktop\NiFeCr_based pure Ni\FeNi\test_nn_effect;
count =0;
% data_result=zeros(12,num_samples);% the initial size of data
fid = fopen(file_path);       
         while ~feof(fid)
          tline = fgetl(fid);
          count= count +1;
          data_result(:,count)= sscanf(tline,'%f'); 
         end
fclose(fid);
read_data=data_result;
end
