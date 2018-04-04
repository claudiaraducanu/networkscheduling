clear all
clear vars

filename     = 'Assignment2.xlsx';

%[~,~,raw] = 
[time] = xlsread(filename,1,'D2:D233');
hour(datetime(time,'ConvertFrom','excel'))
minute(datetime(time,'ConvertFrom','excel'))