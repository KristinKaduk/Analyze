function [dailyStats] = Convert_MatlabOutput_toRInput(monkeyName,dataDir,starting_folder,ending_folder)
% Kristin
% last change 20181106
% Convert_MatlabOutput_toRInput_11062018('Norman','Y:\Data','20180711','20180711')
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
                fname = ['Tst' num2str(year) '-' num2str(months,'%02d') '-' num2str(day,'%02d') '_' num2str(r,'%02d') '.mat'];
                
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
                        x0=[data(i_Trials).task.hnd.tar(i_NrTar).x];
                        y0=[data(i_Trials).task.hnd.tar(i_NrTar).y];
                        r = [data(i_Trials).task.hnd.tar(i_NrTar).radius];
                        
                        is_within = 0;
                        
                        if sqrt(((x0 - x))^2 + (y0 - y)^2) < r % check if (x,y) is within radius r centered on (x0,y0)
                            is_within = 1;
                            %In trials X was the match on the left and the particpant chose left or right?
                            if  x0 < 0 %negativ values are on the left side
                                PosTar_ChoosenOption(i_Trials) = 1;
                            elseif x0 > 0
                                PosTar_ChoosenOption(i_Trials) = 2;
                            end
                        end
                        
                    end
                end
                
                
                
                % time between targ_acq and targ_Hold
                FreqMot_TarHol_to_Sound = repmat(9999,length(nTrials),1); FreqMot_Sound_to_RewardDelivery =repmat(9999,length(nTrials),1); FreqMot_ITI =repmat(9999,length(nTrials),1);
                FreqMot_FixAcq_to_TarHol= repmat(9999,length(nTrials),1); Wagering_ResponseTime= repmat(9999,length(nTrials),1); M2S_ResponseTime = repmat(9999,length(nTrials),1); Freq_RewardDelivery = repmat(9999,length(nTrials),1);
                FreqMot_TarHol_BeforeRewardDelivery= repmat(9999,length(nTrials),1);
                Names= fieldnames(data);
                if strcmp(Names{end-2},'jaw' )
                    [data.('jaw_R')] = data.('jaw');
                    data = rmfield(data,'jaw');
                end
                
                %% plot each trial
                for indTrial = 1: nTrials(end)
                    FreqMot_FixAcq_to_TarHol(indTrial) = NaN;  M2S_ResponseTime(indTrial)  = NaN;Wagering_ResponseTime(indTrial)  = NaN;
                    FreqMot_TarHol_to_Sound(indTrial) = NaN;
                    FreqMot_Sound_to_RewardDelivery(indTrial) = NaN;
                    FreqMot_ITI(indTrial)  = NaN;
                    Freq_RewardDelivery(indTrial) = NaN;
                    FreqMot_TarHol_BeforeRewardDelivery(indTrial) = NaN;
                    %differentiate between correct and incorrect trials
                    if data(indTrial).completed == 1 && data(indTrial).success == 0 %complete & error trial
                        M2S_ResponseTime(indTrial)     =[data(indTrial).states_onset(data(indTrial).states == 5) - data(indTrial).states_onset(data(indTrial).states == 4)];
                         if Task_type(1) == 10
                             Wagering_ResponseTime(indTrial)     =[data(indTrial).states_onset(data(indTrial).states == 24) - data(indTrial).states_onset(data(indTrial).states == 23)];
                             
                         end
                        
                        %1.Motion in the trial fix-acq (state 5) to beginn ITI (state 50)
                        Mot_FixAcq_to_TarHol = data(indTrial).jaw_R(find(data(indTrial).state == 1,1):find(data(indTrial).state == 5,1));
                        FreqMot_FixAcq_to_TarHol(indTrial) = sum(Mot_FixAcq_to_TarHol == 0)/ length(Mot_FixAcq_to_TarHol);
                        %2.2.tar_hol (state 5) to reward (state 21)
                        % the error tone is displayed in state 19 (ABORT) which follows
                        % after the target hold state / the same is when the ITI-Interval
                        % starts
                        Dur_ErrorSound = fid.task.timing.tar_time_hold ;
                        Ind_SoundDisplayed = find(data(indTrial).tSample_from_time_start > (data(indTrial).tSample_from_time_start(find(data(indTrial).state == 5,1))+Dur_ErrorSound),1);
                        Mot_TarHol_to_Sound = data(indTrial).jaw_R(find(data(indTrial).state == 5,1):Ind_SoundDisplayed);
                        FreqMot_TarHol_to_Sound(indTrial) = sum(Mot_TarHol_to_Sound == 0)/ length(Mot_TarHol_to_Sound);
                        %3.reward (state 21) to reward delivery  (state 50)
                        Dur_StartRewardDelivery = fid.task.timing.wait_for_reward + fid.task.timing.tar_time_hold ;
                        Ind_StartRewardDelivery = find(data(indTrial).tSample_from_time_start >= (data(indTrial).tSample_from_time_start(find(data(indTrial).state == 5,1))+Dur_StartRewardDelivery),1);
                        Mot_Sound_to_RewardDelivery =data(indTrial).jaw_R(Ind_SoundDisplayed:Ind_StartRewardDelivery);
                        FreqMot_Sound_to_RewardDelivery(indTrial) = sum(Mot_Sound_to_RewardDelivery == 0)/ length(Mot_Sound_to_RewardDelivery);
                        %4.reward delivery
                        Dur_RewardDelivery = fid.task.reward.time_neutral(1) ;
                        Ind_RewardDelivery = find(data(indTrial).tSample_from_time_start  >= (data(indTrial).tSample_from_time_start(Ind_StartRewardDelivery)+Dur_RewardDelivery),1);
                        Mot_RewardDelivery = data(indTrial).jaw_R(Ind_StartRewardDelivery:Ind_RewardDelivery);
                        Freq_RewardDelivery(indTrial) = sum(Mot_RewardDelivery == 0)/ length(Mot_RewardDelivery);
                        %5.
                        Mot_TarHol_to_BeforeRewardDelivery = data(indTrial).jaw_R(find(data(indTrial).state == 5,1):Ind_StartRewardDelivery-1);
                        FreqMot_TarHol_BeforeRewardDelivery(indTrial) = sum(Mot_TarHol_to_BeforeRewardDelivery == 0)/ length(Mot_TarHol_to_BeforeRewardDelivery);
                        
                        %5.after reward delivery
                        Mot_ITI = data(indTrial).jaw_R(find(data(indTrial).state == 50,1):end);
                        FreqMot_ITI(indTrial) = sum(Mot_ITI == 0)/ length(Mot_ITI);
                        %% Plot the motion data
                        if SavePlot == 1
                            h = figure(indTrial);
                            plot([data(indTrial).tSample_from_time_start],[data(indTrial).jaw_R]);
                            ylim([-0.1, 1.1])
                            line([data(indTrial).tSample_from_time_start(find(data(indTrial).state == 5,1)) data(indTrial).tSample_from_time_start(find(data(indTrial).state == 5,1))], get(gca,'YLim'),'Color','black','LineStyle','--');
                            line([data(indTrial).tSample_from_time_start(find(data(indTrial).state == 50,1)) data(indTrial).tSample_from_time_start(find(data(indTrial).state == 50,1))], get(gca,'YLim'),'Color','black','LineStyle','--');
                            
                            txt1 = 'Start TargHold';
                            text(typecast(double(data(indTrial).tSample_from_time_start(find(data(indTrial).state == 5,1))), 'double'),0.1,txt1)
                            txt1 = 'Start ITI';
                            text(typecast(double(data(indTrial).tSample_from_time_start(find(data(indTrial).state == 50,1))), 'double') + fid.task.timing.wait_for_reward ,0.1,txt1)
                            title(['Norman   ', 'NrTrial:',num2str(indTrial), '   Correct: ',num2str(data(indTrial).success)])
                            pos = get(h,'Position');
                            set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
                            print(h,['Cor_MotionDetection_20180523_Run1_NrTrial',num2str(indTrial) ,'Correct',num2str(data(indTrial).success)],'-dpng','-r0') %dpng
                        end
                        
                    elseif data(indTrial).completed == 1 && data(indTrial).success == 1
                        M2S_ResponseTime(indTrial)     =[data(indTrial).states_onset(data(indTrial).states == 5) - data(indTrial).states_onset(data(indTrial).states == 4)];
                        %M2S_ReleaseTime(indTrial)      =[data(indTrial).states_onset(data(indTrial).states == 5) - data(indTrial).states_onset(data(indTrial).states == 4)];
                         if Task_type(1) == 10
                             Wagering_ResponseTime(indTrial)     =[data(indTrial).states_onset(data(indTrial).states == 24) - data(indTrial).states_onset(data(indTrial).states == 23)];
                         end
                        %1.Motion in the trial fix-acq (state 5) to beginn ITI (state 50)
                        Mot_FixAcq_to_TarHol = data(indTrial).jaw_R(find(data(indTrial).state == 1,1):find(data(indTrial).state == 5,1));
                        FreqMot_FixAcq_to_TarHol(indTrial) = sum(Mot_FixAcq_to_TarHol == 0)/ length(Mot_FixAcq_to_TarHol);
                        %2.tar_hol (state 5) to reward (state 21)
                        %Ind_SoundDisplayed = find(data(insdTrial).tSample_from_time_start == (data(indTrial).tSample_from_time_start(find(data(indTrial).state == 5,1))+ fid.task.timing.tar_time_hold));
                        Mot_TarHol_to_Sound               = data(indTrial).jaw_R(find(data(indTrial).state == 5,1):find(data(indTrial).state == 21,1));
                        FreqMot_TarHol_to_Sound(indTrial) = sum(Mot_TarHol_to_Sound == 0)/ length(Mot_TarHol_to_Sound);
                        %3.reward (state 21) to reward delivery  (state 50)
                        Dur_StartRewardDelivery = fid.task.timing.wait_for_reward ;
                        Ind_StartRewardDelivery = find(data(indTrial).tSample_from_time_start > (data(indTrial).tSample_from_time_start(find(data(indTrial).state == 21,1))+Dur_StartRewardDelivery),1);
                        Mot_Sound_to_StartRewardDelivery          = data(indTrial).jaw_R(find(data(indTrial).state == 21,1):Ind_StartRewardDelivery);
                        FreqMot_Sound_to_RewardDelivery(indTrial) = sum(Mot_Sound_to_StartRewardDelivery == 0)/ length(Mot_Sound_to_StartRewardDelivery);
                        %4.reward delivery
                        Dur_RewardDelivery = fid.task.reward.time_neutral(1) ;
                        Ind_RewardDelivery = find(data(indTrial).tSample_from_time_start > (data(indTrial).tSample_from_time_start(Ind_StartRewardDelivery)+Dur_RewardDelivery),1);
                        Mot_RewardDelivery               = data(indTrial).jaw_R(Ind_StartRewardDelivery:Ind_RewardDelivery);
                        Freq_RewardDelivery(indTrial) = sum(Mot_RewardDelivery == 0)/ length(Mot_RewardDelivery);
                        %5.
                        Mot_TarHol_to_BeforeRewardDelivery = data(indTrial).jaw_R(find(data(indTrial).state == 5,1):Ind_StartRewardDelivery-1);
                        FreqMot_TarHol_BeforeRewardDelivery(indTrial) = sum(Mot_TarHol_to_BeforeRewardDelivery == 0)/ length(Mot_TarHol_to_BeforeRewardDelivery);
                       %5.after reward delivery
                        Mot_ITI               = data(indTrial).jaw_R(find(data(indTrial).state == 50,1):end);
                        FreqMot_ITI(indTrial) = sum(Mot_ITI == 0)/ length(Mot_ITI);
                        
                        %% Plot the motion data
                        if SavePlot == 1
                            h = figure(indTrial);
                            HistPerc  =  hist2per(FreqMot_Sound_to_RewardDelivery);
                            
                            
                            
                            h = figure(indTrial);
                            plot([data(indTrial).tSample_from_time_start],[data(indTrial).jaw_R]);
                            ylim([-0.1, 1.1])
                            line([data(indTrial).tSample_from_time_start(find(data(indTrial).state == 5,1)) data(indTrial).tSample_from_time_start(find(data(indTrial).state == 5,1))], get(gca,'YLim'),'Color','black','LineStyle','--');
                            line([data(indTrial).tSample_from_time_start(find(data(indTrial).state == 21,1)) data(indTrial).tSample_from_time_start(find(data(indTrial).state == 21,1))], get(gca,'YLim'),'Color','black','LineStyle','--');
                            line([data(indTrial).tSample_from_time_start(find(data(indTrial).state == 21,1)) + fid.task.timing.wait_for_reward data(indTrial).tSample_from_time_start(find(data(indTrial).state == 21,1))+ fid.task.timing.wait_for_reward], get(gca,'YLim'),'Color','black','LineStyle','--');
                            
                            txt1 = 'Start TargHold';
                            text(typecast(double(data(indTrial).tSample_from_time_start(find(data(indTrial).state == 5,1))), 'double'),0.1,txt1)
                            txt1 = 'Start Reward';
                            text(typecast(double(data(indTrial).tSample_from_time_start(find(data(indTrial).state == 21,1))), 'double'),0.1,txt1)
                            txt1 = 'Start RewardDelivered';
                            text(typecast(double(data(indTrial).tSample_from_time_start(find(data(indTrial).state == 21,1))), 'double') + fid.task.timing.wait_for_reward ,0.1,txt1)
                            title(['Norman   ', 'NrTrial:',num2str(indTrial), '   Correct: ',num2str(data(indTrial).success)])
                            pos = get(h,'Position');
                            set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
                            print(h,['Cor_MotionDetection_20180523_Run1_NrTrial',num2str(indTrial) ,'Correct',num2str(data(indTrial).success)],'-dpng','-r0') %dpng
                        end %SavePlot
                    end %Completed & Successful
                end %NrTrials
                HistMot_FixAcq_to_ITI = hist2per(FreqMot_FixAcq_to_TarHol);
                HistMot_TarHol_to_Sound = hist2per(FreqMot_TarHol_to_Sound);
                HistMot_Sound_to_RewardDelivery = hist2per(FreqMot_Sound_to_RewardDelivery);
                HistMot_ITI                     = hist2per(FreqMot_ITI);
                
                
                Rotation_Sample = cell2mat(arrayfun(@(data) vertcat(data.hnd.cue(1).shape.rotation),data,'uni',0));
                Rotation_t1 = cell2mat(arrayfun(@(data) vertcat(data.hnd.tar(1).shape.rotation),data,'uni',0));
                Rotation_t2 = cell2mat(arrayfun(@(data) vertcat(data.hnd.tar(2).shape.rotation),data,'uni',0));
                Rot_Diff = abs(Rotation_t1 - Rotation_t2);
                target_selected_Hand = cell2mat(arrayfun(@(data) vertcat(data.target_selected(2)),data,'uni',0));
                target_selected_Eye = cell2mat(arrayfun(@(data) vertcat(data.target_selected(1)),data,'uni',0));%eye, hand ... left = 1< right =2
                
                if Task_type(1) == 10
                    target2_selected_Hand= cell2mat(arrayfun(@(data) vertcat(data.target2_selected(2)),data,'uni',0));
                    target2_selected_Eye= cell2mat(arrayfun(@(data) vertcat(data.target2_selected(1)),data,'uni',0));%e
                end
                %% ABORT CODE
                abort_code = arrayfun(@(data) vertcat(data.abort_code),data,'uni',0);
                T = [];
                for indTrials =  1: length(data)
                    T.abort_code(1,indTrials) = {abort_code{indTrials}}';
                end
                T.abort_code =        T.abort_code' ;
                Abort_code = struct2table(T);
                %%
                 switch  Task_type(1)
                    case 10
                        
                        Tab_Trial = array2table([nTrials; completed; success ; target_selected_Hand; target_selected_Eye; PosTar_ChoosenOption; target2_selected_Hand;...
                            target2_selected_Eye; aborted_state; Task_type;Rotation_Sample; Rotation_t1; Rotation_t2; Rot_Diff;...
                            reward_selected; reward_time; rewarded; FreqMot_FixAcq_to_TarHol';HistMot_FixAcq_to_ITI'; FreqMot_TarHol_to_Sound';HistMot_TarHol_to_Sound'; FreqMot_Sound_to_RewardDelivery';...
                            Freq_RewardDelivery'; FreqMot_ITI';FreqMot_TarHol_BeforeRewardDelivery'; HistMot_ITI';M2S_ResponseTime'; Wagering_ResponseTime']');
                        % give the colums a name
                        Tab_Trial.Properties.VariableNames = {getname(nTrials),getname(completed),getname(success),getname(target_selected_Hand),...
                            getname(target_selected_Eye),  getname(PosTar_ChoosenOption),getname(target2_selected_Hand), getname(target2_selected_Eye),getname(aborted_state) ,getname(Task_type),getname(Rotation_Sample),getname(Rotation_t1),...
                            getname(Rotation_t2),getname(Rot_Diff),getname(reward_selected),getname(reward_time),getname(rewarded),getname(FreqMot_FixAcq_to_TarHol),getname(HistMot_FixAcq_to_ITI),...
                            getname(FreqMot_TarHol_to_Sound), getname(HistMot_TarHol_to_Sound),getname(FreqMot_Sound_to_RewardDelivery),...
                            getname(Freq_RewardDelivery),getname(FreqMot_ITI),getname(FreqMot_TarHol_BeforeRewardDelivery),getname(HistMot_ITI),getname(M2S_ResponseTime),getname(Wagering_ResponseTime)};
                        
                    case 9
                        % create a table
                        Tab_Trial = array2table([nTrials; completed; success ; target_selected_Hand; target_selected_Eye; PosTar_ChoosenOption; aborted_state; Task_type;Rotation_Sample; Rotation_t1; Rotation_t2; Rot_Diff;...
                            FreqMot_FixAcq_to_TarHol';HistMot_FixAcq_to_ITI'; FreqMot_TarHol_to_Sound';HistMot_TarHol_to_Sound'; FreqMot_Sound_to_RewardDelivery';...
                            Freq_RewardDelivery'; FreqMot_ITI';FreqMot_TarHol_BeforeRewardDelivery';HistMot_ITI';M2S_ResponseTime']');
                        % give the colums a name
                        Tab_Trial.Properties.VariableNames = {getname(nTrials),getname(completed),getname(success),getname(target_selected_Hand),getname(target_selected_Eye),...
                            getname(PosTar_ChoosenOption), getname(aborted_state) ,getname(Task_type),getname(Rotation_Sample),getname(Rotation_t1),...
                            getname(Rotation_t2),getname(Rot_Diff),getname(FreqMot_FixAcq_to_TarHol),getname(HistMot_FixAcq_to_ITI),...
                            getname(FreqMot_TarHol_to_Sound), getname(HistMot_TarHol_to_Sound),getname(FreqMot_Sound_to_RewardDelivery),...
                            getname(Freq_RewardDelivery),getname(FreqMot_ITI),getname(FreqMot_TarHol_BeforeRewardDelivery),getname(HistMot_ITI),getname(M2S_ResponseTime)};
                 end
                
                
                
                DataOneFile  = [];
                T  = [];
                Monkey  = fname(1:3);
                Date    = fname(4:13);
                Run     = fname(15:16);
                T = table(repmat(Monkey,size(Tab_Trial(:,1),1),1), repmat(Date,size(Tab_Trial(:,1),1),1),repmat(Run,size(Tab_Trial(:,1),1),1) ,'VariableNames', {'Monkey' 'Date' 'Run'} );
                DataOneFile = [T, Tab_Trial, Abort_code];
                
                Data = [ DataOneFile;Data];
                
                
            elseif Task_type(1) ==2 &&  strcmp(monkeyName,'Cornelius')
                
                nTrials = []; reach_hand = [];Tab_Trial= [];
                nTrials       =[data.n];
                for indTrial = 1: nTrials(end)
                    if  isempty(data(indTrial).reach_hand)
                        reach_hand(indTrial) = NaN;
                    else
                        reach_hand(indTrial) =   [data(indTrial).reach_hand(1)];
                    end
                end
                force_hand    =[data.force_hand];
                completed     =[data.completed];
                effector      =[data.effector];
                rest_hand     =[data.rest_hand];
                hand_choice   =[data.hand_choice];
                choice        =[data.choice];
                aborted_state =[data.aborted_state];
                success       =[data.success];
                Task_type     =[data.type];
                target_selected_Hand = cell2mat(arrayfun(@(data) vertcat(data.target_selected(2)),data,'uni',0));
                target_selected_Eye = cell2mat(arrayfun(@(data) vertcat(data.target_selected(1)),data,'uni',0));
                
                Tab_Trial = array2table([nTrials;          completed;         success; reach_hand; target_selected_Hand; target_selected_Eye; aborted_state; Task_type]'  ,...
                    'VariableNames',{getname(nTrials),  getname(completed),getname(success),getname(reach_hand),getname(target_selected_Hand),getname(target_selected_Eye),getname(aborted_state) ,getname(Task_type)});
                
                DataOneFile  = [];
                T  = [];
                Monkey  = fname(1:3);
                Date    = fname(4:13);
                Run     = fname(15:16);
                T = table(repmat(Monkey,size(Tab_Trial(:,1),1),1), repmat(Date,size(Tab_Trial(:,1),1),1),repmat(Run,size(Tab_Trial(:,1),1),1) ,'VariableNames', {'Monkey' 'Date' 'Run'} );
                %DataOneFile = [T, Tab_Trial, ErrornameTableColum];
                DataOneFile = [T, Tab_Trial];
                Data = [ DataOneFile;Data];
                
            end %TaskType
            
            
            
            
            
        end
    end
    
    
end
writetable(Data,['Y:\Projects\Wagering_monkey\Data\VariablesInTableForR\' monkeyName, '\',Monkey,  'M2S_psychophysicTask_since' ,starting_folder,'_until_', ending_folder,'.txt'], 'Delimiter', ',')
writetable(Data,['C:\Users\kkaduk\Dropbox\promotion\Projects\Wagering_Monkey\Results\' Monkey, '\',Monkey,  'M2S_psychophysicTask_since' ,starting_folder,'_until_', ending_folder,'.txt'], 'Delimiter', ',')
disp(['saved  C:\Users\kkaduk\Dropbox\promotion\Projects\Wagering_Monkey\Results\' monkeyName, '\',Monkey,  'M2S_psychophysicTask_since' ,starting_folder,'_until_', ending_folder,'.txt'])
end


