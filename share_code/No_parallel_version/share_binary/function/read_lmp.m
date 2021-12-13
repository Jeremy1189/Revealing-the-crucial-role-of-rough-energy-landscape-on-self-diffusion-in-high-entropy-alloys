function [data,infor_set] = read_lmp(lmp_file_path)
data=[];
count =0;
fid = fopen(lmp_file_path);       
        flag_char = {'Atoms # atomic'};
        while ~feof(fid)
          tline = fgetl(fid);
          count= count +1;
          infor_set{count} = tline; %#ok<AGROW>
          for num_char = 1: length(flag_char )
           matches = strfind(tline, flag_char{num_char});
           num = length(matches);
           if num > 0 && matches == 1
               tline = fgetl(fid);%empty line
               while ~feof(fid)
                 tline = fgetl(fid);
                 data = [data,sscanf(tline, '%f %f %f %f %f')]; %#ok<AGROW>   
               end
           end           
          end
%           tline = fgetl(fid);
        end
    fclose(fid);
end