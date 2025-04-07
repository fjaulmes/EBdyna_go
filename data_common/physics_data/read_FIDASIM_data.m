filename='neut_rate_table.bin'
fileID = fopen(filename,'w');
fwrite(fileID,[1:9]);
fclose(fileID);

%%
fileID = fopen(filename, 'r', 'ieee-le');
if fileID == -1, error('Cannot open file: %s', filename); end
format = 'uint16';
Data = fread(fileID, Inf, format);
fclose(fileID);