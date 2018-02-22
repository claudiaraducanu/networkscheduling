%%  Initialization
%addpath('D:Documents\cplex\matlab\x64_win64'); %Luka
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1271\cplex\matlab\x64_win64'); %Satrio
%addpath('Users\Luka\Documents\IBM\ILOG\CPLEX_Studio1271'); %Bryan
clc
clearvars
close all
warning('off','MATLAB:lang:badlyScopedReturnValue')
warning('off','MATLAB:xlswrite:NoCOMServer')


%%  Determine input
%   Select input file and sheet
input        =   [pwd '/Input_AE4424_Ass1P1.xlsx'];
%filsol      =   'Solutions.xlsx';

cost = xlsread(input,1,'D2:D31');
capa = xlsread(input,1,'E2:E31');
cost = xlsread(input,2,'D2:D41');



