%% Initialize
clearvars
clear all

%% Read file
filename     = 'Input_AE4424_Ass1P1.xlsx';
[~, dname]   = xlsread(filename,1,'B2:B31'); % names of departure airports
[~, aname]   = xlsread(filename,1,'C2:C31'); % names of arrival airports
acost        = xlsread(filename,1,'D2:D31'); % cost associated to existing arc
acap         = xlsread(filename,1,'E2:E31'); % cost associated to existing arc
 
edges       = [ dname, aname];

nodenumber     = 1;               % first node equal to first origin node;
orgname        = edges{1,1};      % name of first node;
edges{1,3}     = nodenumber;      % set number of first node;


for i = 1:(size(edges,1))
    if  strcmp(edges{i,1},orgname) == 0
        nodenumber    = nodenumber+1;
        edges{i,3} = nodenumber;
        orgname   = edges{i,1};
    else
        edges{i,3} = nodenumber;
    end
end

for i = 1:(size(edges,1))
   if  sum(strcmp(edges{i,2},edges(:,1))) > 0
       indx = find(strcmp(edges{i,2},edges(:,1)));
       edges{i,4} = edges{indx(1),3};
   else
       indx = find(strcmp(edges{i,2});
       nodenumber = nodenumber+1;
       edges{i,4} = nodenumber;
   end
end

