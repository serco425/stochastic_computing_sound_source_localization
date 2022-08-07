% This file takes some simulated audio data and creates a c header file to
% be included for testing purposes.

addpath('./stochastic_computing/');
load('simulation_data.mat','X','stdvalues');    % get some audio data
software_src_path='./../LA_C_Code/C_Cross_Correlation/'; % path where testSoftware is Located

X_scaled = X.*50;  %do some scaling
N_print = length(X_scaled);

%check if scaling results in xcorr lower than 1;
% figure();plot(xcorr(X_scaled(:,1),X_scaled(:,2)));


bit_length = 8;

type = numerictype(1,bit_length,bit_length-1);

fi_left = fi(X_scaled(:,1),type);
fi_right = fi(X_scaled(:,2),type);
% figure();plot(xcorr(fi_left,fi_left,3));

 

fid = fopen([software_src_path 'AUDIO_data.h'],'w');

fprintf(fid,'// Include File to Test Cross-Correlation-Software\r\n');
fprintf(fid,'//created on %s\r\n',date);
fprintf(fid,'//by file \LA_matlab\audio_to_c_file.m\r\n');
fprintf(fid,'//represents audio data of one microphone\r\n');
fprintf(fid,'//- sample rate is 12kHz\r\n');
fprintf(fid,'#ifndef AUDIO_DATA_H\r\n');
fprintf(fid,'#define AUDIO_DATA_H\r\n');

fprintf(fid,sprintf("uint8_t AUDIO_data_left[%d] = {\r\n",N_print));
ind = 1;
while (ind <= N_print-1)

fprintf(fid,'%d,',fi_left.int(ind));

%for better readability insert some newlines
if (mod(ind,100)==0)
    fprintf(fid,'\r\n');
end

ind = ind+1;
end
fprintf(fid,'%d};\r\n',fi_left.int(end));

fprintf(fid,sprintf("uint8_t AUDIO_data_right[%d] = {\r\n",N_print));
ind = 1;
while (ind <= N_print-1)

fprintf(fid,'%d,',fi_right.int(ind));

%for better readability insert some newlines
if (mod(ind,100)==0)
    fprintf(fid,'\r\n');
end

ind = ind+1;
end
fprintf(fid,'%d};\r\n',fi_right.int(end));

fprintf(fid,'#endif\r\n');
fclose(fid);