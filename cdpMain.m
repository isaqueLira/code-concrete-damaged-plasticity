%% Opening
clear % clean the workspace
clc % clean the command window

%% Input Data
inputData = inputdlg({'fck (MPa): ','Length MEF (mm): ','b: '});

fck=str2double(inputData{1}); 
leq=str2double(inputData{2}); 
b=str2double(inputData{3});

%% Results
[RC,RT,DC,DT,arq] = results(fck,leq,b);