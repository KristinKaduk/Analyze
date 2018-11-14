function [dailyStats] = Convert_MatlabOutput_toRInput(monkeyName,dataDir,starting_folder,ending_folder)

 %  Convert_MatlabOutput_toRInput('Norman','D:\Data','20180323','20180323')
 %  Convert_MatlabOutput_toRInput('Norman','Y:\Data','20180323','20180323')
[number]    =number_of_folders(dataDir,starting_folder,ending_folder);
startDate   =starting_folder;
yearMonth   =starting_folder(1:6);
endDate     =ending_folder;
f_name      =num2str(starting_folder);
last_two    =f_name(end-1:end); 
Data        =[]; 
for j=1:number

% if nargin<1
%     monkeyName = 'Tesla';
% end
%   dataDir ='C:\Users\M.Koester\Desktop\LindaNieuw';
 % dataDir =['Y:\Magnus\201305' num2str(last_two) filesep];
%   dataDir =['X:\Magnus\201303' num2str(last_two) filesep]; %Danial
    dataDir =[dataDir ,filesep,monkeyName,filesep, yearMonth num2str(last_two) filesep]; %Danial setup 2




% dataDir = [upper(dataDriveLetter) ':/' monkeyName '/'];   % specify data directory
%dataDir = [dataDir filesep];   % specify data directory
if ~isdir(dataDir); error('No valid data directory specified'); end

if nargin < 5
    y1=2011; m1=01; d1=01;  % start date
    y2=2020; m2=01; d2=01;  % end date
elseif numel(startDate)==5;
    y1=startDate(1); m1=startDate(2); d1=startDate(3);
    y2=endDate(1);   m2=endDate(2);   d2=endDate(3);
else
    startDate=str2num(startDate);    endDate=str2num(endDate); %Danial
    y1=floor(startDate/1e4);            % convert 8 digit number to date
    m1=floor((startDate-y1*1e4)/1e2);
    d1=startDate-y1*1e4-m1*1e2;
    y2=floor(endDate/1e4);
    m2=floor((endDate-y2*1e4)/1e2);
    d2=endDate-y2*1e4-m2*1e2;
end


%% 1) Find all available txt files
% List all availabe .txt files from dataDir
% runsAvailable is an (nMatFiles x [year month day run]) array
D = dir(dataDir);
runsAvailable = nan(length(D),4);   % preallocate
c = 0;
for i = 3:length(D)
    if strcmp(D(i).name(end-2:end),'mat')   % if .mat file found
        c = c+1;
        runsAvailable(c,1)  = str2double(D(i).name(4:7));   % year
        runsAvailable(c,2)  = str2double(D(i).name(9:10));  % month
        runsAvailable(c,3)  = str2double(D(i).name(12:13)); % day
        runsAvailable(c,4)  = str2double(D(i).name(15:16)); % run
    end
end
runsAvailable = runsAvailable(1:c,:); % truncate to actual size

datesAvailable = runsAvailable( logical([1; diff(runsAvailable(:,3))]) ,1:3);   % all dates for which data is availabel (n x [year month day])

% check which dates fall into speficied date range
isAfter =  (datesAvailable(:,1) > y1 | ...
            datesAvailable(:,1) == y1 & datesAvailable(:,2) > m1 | ...
            datesAvailable(:,1) == y1 & datesAvailable(:,2) == m1 & datesAvailable(:,3) >= d1);
isBefore=  (datesAvailable(:,1) < y2 | ...
            datesAvailable(:,1) == y2 & datesAvailable(:,2) < m2 | ...
            datesAvailable(:,1) == y2 & datesAvailable(:,2) == m2 & datesAvailable(:,3) <= d2);        
withinRange = isAfter & isBefore;
datesAvailable(~withinRange,:) = [];

%% 2) Loop over all dates for which data is available. For each date loop over all runs.
% dailyStats: (day x [year month day nTrials hitRate])
nDays = size(datesAvailable,1);
dailyStats = nan(nDays,8);
dailyTimeReward = nan(nDays,1);
dailyWeekDay = repmat(char(0),nDays,3);
for i = 1:nDays
    y = datesAvailable(i,1);
    m = datesAvailable(i,2);
    d = datesAvailable(i,3);
    runs = runsAvailable(runsAvailable(:,1)==y & runsAvailable(:,2)==m & runsAvailable(:,3)==d,4);

    for r = runs';            
        fname = [monkeyName(1:3) num2str(y) '-' num2str(m,'%02d') '-' num2str(d,'%02d') '_' num2str(r,'%02d') '.mat'];

        fid=load([dataDir fname]);
        data=fid.trial;
nTrials = []; reach_hand = []; 
nTrials       =[data.n];
completed     =[data.completed];
aborted_state =[data.aborted_state];
success       =[data.success];
Task_type     =[data.type];
Rotation_t1 = cell2mat(arrayfun(@(data) vertcat(data.hnd.tar(1).shape.rotation),data,'uni',0));
Rotation_t2 = cell2mat(arrayfun(@(data) vertcat(data.hnd.tar(2).shape.rotation),data,'uni',0));
Rot_Diff = abs(Rotation_t1 - Rotation_t2);
target_selected = cell2mat(arrayfun(@(data) vertcat(data.target_selected(2)),data,'uni',0));
abort_code = arrayfun(@(data) vertcat(data.abort_code),data,'uni',0);

getname = @(x) inputname(1);
Tab_Trial = array2table([nTrials; completed; success; aborted_state; Task_type; Rotation_t1; Rotation_t2; Rot_Diff; target_selected]'  ,...
    'VariableNames',{getname(nTrials),getname(completed),getname(success),getname(aborted_state),getname(Task_type),getname(Rotation_t1),...
    getname(Rotation_t2),getname(Rot_Diff),getname(target_selected)});

DataOneFile  = [];
T  = [];
Monkey  = fname(1:3);
Date    = fname(4:13);
Run     = fname(15:16);     
T = table(repmat(Monkey,size(Tab_Trial(:,1),1),1), repmat(Date,size(Tab_Trial(:,1),1),1),repmat(Run,size(Tab_Trial(:,1),1),1) ,'VariableNames', {'Monkey' 'Date' 'Run'} );
DataOneFile = [T, Tab_Trial];
Data = [ DataOneFile;Data];
    end       
    end

    
end
writetable(Data,['Y:\Projects\Wagering_monkey\Results\Norman\',Monkey,  'M2S_psychophysicTask_since' ,starting_folder,'_until_', ending_folder,'.txt'], 'Delimiter', ' ')
writetable(Data,['C:\Users\kkaduk\Dropbox\promotion\Projects\Wagering_Monkey\Results\Norman\',Monkey,  'M2S_psychophysicTask_since' ,starting_folder,'_until_', ending_folder,'.txt'], 'Delimiter', ' ')
end


