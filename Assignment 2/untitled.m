clear all
clear vars

filename     = 'Assignment2.xlsx';

% import origin and destination airport names

[~,org]  = xlsread(filename,1,'B2:B233'); % origin
[~,dest] = xlsread(filename,1,'C2:C233'); % destimations

L        = size(org,1); % number of flight arcs

% locations in time-space network

[a,~,X]  = unique(org); % airport names
A        = size(a,1);   % number of airports 

% represent times in a day in minutes from 0 (0:00) to 1439 (23:59)  

dep_time     = xlsread(filename,1,'D2:D233');
[h,m]        = hms(datetime(dep_time,'ConvertFrom','excel'));
dep_time     = h*60+m;            % departure time in minutes

arr_time     = xlsread(filename,1,'E2:E233'); 
[h,m]        = hms(datetime(arr_time,'ConvertFrom','excel'));
arr_time     = h*60+m;            % arrival time in minutes

% aircraft data per type
[~,varnames] = xlsread(filename,4,'A1:D1');
[~,~,AC]   = xlsread(filename,4,'A2:D5');
AC = cell2table(AC,'VariableNames',varnames);

% cost of operating aircraft of type k on fight i 
cost         = xlsread(filename,1,'F2:I233');

% logical array, where 1 represents the fact that that flight
% can be operated by aircraft type k
idx_fl       = ~isnan(cost); 

% Flight arcs for each aircraft. 

FL.A330 = cell2table(sortrows([org num2cell(dep_time) dest num2cell(arr_time + AC.TAT(1))],[1,2]),...
            'VariableNames',{'Org' 'Dept' 'Dest' 'Arr'});
FL.A340 = cell2table(sortrows([org num2cell(dep_time) dest num2cell(arr_time + AC.TAT(2))],[1,2]),...
            'VariableNames',{'Org' 'Dept' 'Dest' 'Arr'});
FL.B737 = cell2table(sortrows([org num2cell(dep_time) dest num2cell(arr_time + AC.TAT(3))],[1,2]),...
            'VariableNames',{'Org' 'Dept' 'Dest' 'Arr'});
FL.A738 = cell2table(sortrows([org num2cell(dep_time) dest num2cell(arr_time + AC.TAT(4))],[1,2]),...
            'VariableNames',{'Org' 'Dept' 'Dest' 'Arr'});

% Ground arc determination

arc     = zeros(A,1439);

Y = hist(X,unique(X));
for j = 1:A
   idx          = find(strcmp(FL.A330.Org, a(j))); % index of nodes departing at airport j
   arc(j,FL.A330.Dept(idx)') = 1;                  % if node is departing set 1 
   idx          = find(strcmp(FL.A330.Org, a(j))); % index of nodes departing at airport j
   arc(j,FL.A330.Dept(idx)') = 1;                  % if node is departing set 1 
end

%isnan(cost)


%cost = table(xlsread(filename,1,'F2:F233'),xlsread(filename,1,'G2:G233'),...
%    xlsread(filename,1,'H2:H233'),xlsread(filename,1,'I2:I233'),...
%    'VariableNames',varnames);


% flight arc determination for each aircraft

% create table of flight arcs. 

% nodes.depart  = sortrows([ org num2cell(dep_time_min) ],[1,2]); % node 
% nodes.arr     = sortrows([ dest num2cell(arr_time_min)],[1,2]); % node
% 
% % aircraft data per type




