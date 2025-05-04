%% getting from KE conserving case
clear
% close all
% filename = "simulation_default_rerun";
filename = "timing/timing_de";
startRow = 3;
formatSpec = '%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
simulation = [dataArray{1:end-1}];
clearvars filename startRow formatSpec fileID dataArray ans;
start_index=10;
time_de=simulation(start_index:end-1,2);
totaltime_de=simulation(start_index:end-1,3);
VOF_de=simulation(start_index:end-1,4);
Velpred_de=simulation(start_index:end-1,5);
Press_de=simulation(start_index:end-1,6);

% filename = "simulation_ke_case4_new";
filename = "timing/timing_ke";
startRow = 3;
formatSpec = '%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
simulation = [dataArray{1:end-1}];
clearvars filename startRow formatSpec fileID dataArray ans;
start_index=10;
time_KE=simulation(start_index:end-1,2);
totaltime_KE=simulation(start_index:end-1,3);
VOF_KE=simulation(start_index:end-1,4);
Velpred_KE=simulation(start_index:end-1,5);
Press_KE=simulation(start_index:end-1,6);



%%
fprintf('De Total   time average: %.16f \n',mean(totaltime_de))

fprintf('De VOF     time average: %.16f \n',mean(VOF_de))

fprintf('De VelPred time average: %.16f \n',mean(Velpred_de))

fprintf('De Press   time average: %.16f \n',mean(Press_de))


fprintf('KE Total   time average: %.16f \n',mean(totaltime_KE))

fprintf('KE VOF     time average: %.16f \n',mean(VOF_KE))

fprintf('KE VelPred time average: %.16f \n',mean(Velpred_KE))

fprintf('KE Press   time average: %.16f \n',mean(Press_KE))




