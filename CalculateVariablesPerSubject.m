function [dailyStats] = CalculateVariablesPerSubject(monkeyName,dataDir,starting_folder,ending_folder)
% Kristin
% last change 20181106
% CalculateVariablesPerSubject('Norman','Y:\Data','20181130','20181130')
%% TODO
%remove dates
%% PARAMETER
SavePlot = 0;



%%
startDate   =starting_folder;
yearMonth   =starting_folder(1:6);
endDate     =ending_folder;
f_name      =num2str(starting_folder);
last_two    =f_name(end-1:end);
Data        =[];
indFile  =[];
Direct =[dataDir ,filesep,monkeyName,filesep];
Files = dir(Direct);  %% list of all folders
for i = 1:length(Files)
    if strcmp(Files(i).name, starting_folder)==1 ||  strcmp(Files(i).name, ending_folder)==1
        indFile = [indFile, i];
    end
end
% delete all Folders which are not used
if size(indFile,2) == 1
    Files = Files(indFile(1));
else
    Files = Files(indFile(1):indFile(2));
end
[number] = size(Files,1);

for j=1:number
    
    %% assumption that the datafiles are contiously in number (What happends with weekends?)
    %     if j == 1
    %    Dir =[Direct Files(indFile(j)).name filesep];
    %     elseif j == number
    %    Dir =[Direct Files(indFile(2)).name filesep];
    %     else
    %     Dir =[Direct num2str(str2double(Files(indFile(1)).name)+j-1) filesep];
    %     end
    
    Dir =[Direct Files(j).name filesep];
    
    
    if ~isdir(Dir); error('No valid data directory specified'); end
    
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
    D = dir(Dir);
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
        year = datesAvailable(i,1);
        months = datesAvailable(i,2);
        day = datesAvailable(i,3);
        runs = runsAvailable(runsAvailable(:,1)==year & runsAvailable(:,2)==months & runsAvailable(:,3)==day,4);
        
        for r = runs';
            
            fname = [monkeyName(1:3) num2str(year) '-' num2str(months,'%02d') '-' num2str(day,'%02d') '_' num2str(r,'%02d') '.mat'];
            
            if strcmp(monkeyName, 'Test')
                fname = ['Tst' num2str(y) '-' num2str(months,'%02d') '-' num2str(day,'%02d') '_' num2str(r,'%02d') '.mat'];
                
            end
            
            fid=load([Dir fname]);
            data=fid.trial;
            Task_type     =[data.type];
            getname = @(x) inputname(1);
            
            if Task_type(1) == 10  ||Task_type(1) == 9
                
                
                nTrials = []; reach_hand = [];aborted_state = [];completed= [];
                nTrials       =[data.n];
                force_hand    =[data.force_hand];
                completed     =[data.completed];
                effector      =[data.effector];
                rest_hand     =[data.rest_hand];
                reach_hand    =[data.reach_hand];
                hand_choice   =[data.hand_choice];
                choice        =[data.choice];
                aborted_state =[data.aborted_state];
                success       =[data.success];
                Task_type     =[data.type];
                if Task_type(1) == 10
                    reward_selected  = [data.reward_selected ];
                    reward_time = [data.reward_time ];
                    rewarded    = [data.rewarded ];
                end
                  target_selected_Hand= cell2mat(arrayfun(@(data) vertcat(data.target_selected(2)),data,'uni',0));
                    target_selected_Eye= cell2mat(arrayfun(@(data) vertcat(data.target_selected(1)),data,'uni',0));%e
                

                    target2_selected_Hand= cell2mat(arrayfun(@(data) vertcat(data.target2_selected(2)),data,'uni',0));
                    target2_selected_Eye= cell2mat(arrayfun(@(data) vertcat(data.target2_selected(1)),data,'uni',0));%e
               %% read out the position for completed trials
                %  x,z coordinate of touch the screen during the target
                %  acquisition period
                
                PosTar_ChoosenOption = nan(1,length(data));
                for i_Trials = find([data.completed]== 1)
                    x_hnd = [data(i_Trials).x_hnd]';                  y_hnd = [data(i_Trials).y_hnd]';
                    StatesPerTrial = [data(i_Trials).state]';
                    IndStateTarAcq = find(StatesPerTrial == 4); % state == 4 tar_acq
                    x= x_hnd(IndStateTarAcq(end));
                    y= y_hnd(IndStateTarAcq(end));
                    
                    % x,y coordinate of the targets
                    % Which target has been touched? - left or right?
                    for i_NrTar = 1: length(data(i_Trials).task.hnd.tar) %target 1 is the match
                        x0=[data(i_Trials).task.hnd.tar(i_NrTar).x]; %x-position of the target1 etc.
                        y0=[data(i_Trials).task.hnd.tar(i_NrTar).y];
                        r = [data(i_Trials).task.hnd.tar(i_NrTar).radius];
                        
                        is_within = 0;
                        
                        if sqrt(((x0 - x))^2 + (y0 - y)^2) < r % check if (x,y) is within radius r centered on (x0,y0)
                            is_within = 1;
                            %left or right?
                            if  x0 < 0 %negativ values are on the left side
                                PosTar_ChoosenOption(i_Trials) = 1; %match presented on the left 
                            elseif x0 > 0
                                PosTar_ChoosenOption(i_Trials) = 2; %match presented on the right 
                            end
                        end
                        
                    end
                end                
n_wagers = unique(target2_selected_Hand);
n_wagers = n_wagers(~isnan(unique(target2_selected_Hand)));
WageringVSControl = ones(length(completed),1);
%target 1 is the match - was the match left or right
Response_Left = target_selected_Hand;
%Stimulus_Left = ;
% out = metaD_PerSubject(n_wagers, WageringVSControl, completed,selected_Left, target2_selected_Hand,selected_Right);
% General(file_index).d       = out.da;
% General(file_index).metaD   = out.meta_da;
% General(file_index).metaEfficiency = out.M_ratio;
% General(file_index).Typ1_criterion = out.Typ1_criterion;
%% SLOPE-BASED METACOGNITION
% linear fit of all correct trials as function of Wagers of one subject
        p_conf_post      = polyfit(1:n_wagers,PostCertain.PercSuc,1);
		r_conf_post      = p_conf_post(1) .* [1:n_wagers] + p_conf_post(2);
		p_err_post       = polyfit(1:n_wagers,PostCertain.PercUnsuc,1);
		r_err_post       = p_err_post(1) .* [1:n_wagers] + p_err_post(2);
		
        p_conf_pre       = polyfit(1:n_wagers,PreCertain.PercSuc,1);
		r_conf_pre       = p_conf_pre(1) .* [1:n_wagers] + p_conf_pre(2);
		p_err_pre        = polyfit(1:n_wagers,PreCertain.PercUnsuc,1);
		r_err_pre        = p_err_pre(1) .* [1:n_wagers] + p_err_pre(2);
        r_conf_corrected = r_conf_post(:) - r_conf_pre(:);
		r_err_corrected  = r_err_post(:)  - r_err_pre(:);
		
		meta_conf          = p_conf_post(1) - p_conf_pre(1);
		meta_err           = -p_err_post(1) + p_err_pre(1);
		meta_conf_pre      = p_conf_post(1) - p_conf_pre(1);
		meta_err_pre       = -p_err_post(1) + p_err_pre(1);
		General(file_index).meta_slope               = meta_conf(:) + meta_err(:);
		meta_pre           = meta_conf_pre + meta_err_pre;
        
            end
        end
    end
end
writetable(Data,['Y:\Projects\Wagering_monkey\Results\' monkeyName, '\',Monkey,  'Variables_perSubject_since' ,starting_folder,'_until_', ending_folder,'.txt'], 'Delimiter', ',')
writetable(Data,['C:\Users\kkaduk\Dropbox\promotion\Projects\Wagering_Monkey\Results\' Monkey, '\',Monkey,  'Variables_perSubject_since' ,starting_folder,'_until_', ending_folder,'.txt'], 'Delimiter', ',')
disp(['saved  C:\Users\kkaduk\Dropbox\promotion\Projects\Wagering_Monkey\Results\' monkeyName, '\',Monkey,  'Variables_perSubject_since' ,starting_folder,'_until_', ending_folder,'.txt'])
end
