function ComputeVariablesPerSession_M2S_Wagering(path_Data, path_save, dateStart, dateEnd, monkey)

% INPUT
% path_Data = 'Y:\Projects\Wagering_monkey\Data\VariablesInTableForR'; 
% path_save = 'Y:\Projects\Wagering_monkey\Data\VariablesInTableForR'; 
% dateStart = '20190103'; 
% dateEnd = '20190111'; 
% monkey = 'Norman'; 
%%

load([path_Data, filesep, monkey, filesep, monkey(1:3),'M2S_Wagering_since',dateStart,'_until_' dateEnd ])
for i = 1: size(Data,1)
Date(i) = {Data.Date(i,:)}; 
end

DifferentDates = sort(unique(Date)); 
AllData_Combined = []; 
General = []; 
Tab_Session_Diff =[];
for iDate = 1: size(DifferentDates,2)
    % find all entries for a specific Date
 Table = []; 
    for iSubset = 1:  size(Data,1)
    if all(Data.Date(iSubset,:) == char(DifferentDates(iDate))) ==1 
    Table = [Table;Data(iSubset ,:)]; 
    end
    end
General(iDate).Monkey     = monkey(1:3);    
General(iDate).Date       = DifferentDates(iDate);
%% Meta-D
wagers                                  = unique(Table.target2_selected_Hand); 
wagers(isnan(wagers))                 = []; 
n_wagers = max(wagers); 
selected_Right                            =Table.PosTar_ChoosenOption';
selected_Left                             =Table.PosTar_ChoosenOption';
selected_Right(selected_Right == 1) = NaN; 
selected_Right(selected_Right == 2) = 1; 
selected_Left(selected_Left == 2) = NaN; 

    
              
out = metaD_PerSubject(n_wagers, repmat(0,1,length( Table.success)), Table.success',selected_Left, Table.target2_selected_Hand',selected_Right);
General(iDate).d       =  out.da; 
General(iDate).metaD   = out.meta_da;
General(iDate).metaEfficiency = out.M_ratio;
General(iDate).metaEfficiency_Difference = out.M_diff;
General(iDate).meta_c1       =  out.type2_fit.meta_c1; 
General(iDate).c_1            =  out.c_1; 


%% Meta-D per Difficulty Level
TDiff = []; Tab_Diff =[];
out= [];
for iDiff = 1: length(unique(Table.Rot_Diff))
Diff = sort(unique(Table.Rot_Diff)) ;
DiffTable = Table(Table.Rot_Diff ==  Diff(iDiff),:);
selected_Right = [];selected_Left = [];n_wagers =[];
wagers                                  = unique(DiffTable.target2_selected_Hand); 
wagers(isnan(wagers))                 = []; 
n_wagers = max(wagers); 
selected_Right                            =DiffTable.PosTar_ChoosenOption';
selected_Left                             =DiffTable.PosTar_ChoosenOption';
selected_Right(selected_Right == 1) = NaN; 
selected_Right(selected_Right == 2) = 1; 
selected_Left(selected_Left == 2) = NaN; 
out = metaD_PerSubject(n_wagers, repmat(0,1,length( DiffTable.success)), DiffTable.success',selected_Left, DiffTable.target2_selected_Hand',selected_Right);
TDiff(iDiff).Monkey     = monkey(1:3);    
TDiff(iDiff).Date       = DifferentDates(iDate);
TDiff(iDiff).DifficultyLevel= Diff(iDiff);
if out.da ==Inf

TDiff(iDiff).d                =  NaN; 
TDiff(iDiff) .metaD           = NaN;  
TDiff(iDiff).metaEfficiency   = NaN; 
TDiff(iDiff).metaEfficiency_Difference = NaN; 
TDiff(iDiff).meta_c1          =  NaN;  
TDiff(iDiff).c_1              =  NaN;  
else
TDiff(iDiff).d               =  out.da; 
TDiff(iDiff) .metaD           = out.meta_da;
TDiff(iDiff).metaEfficiency   = out.M_ratio;
TDiff(iDiff).metaEfficiency_Difference = out.M_diff;
TDiff(iDiff).meta_c1       =  out.type2_fit.meta_c1; 
TDiff(iDiff).c_1            =  out.c_1; 
end
end

%% SLOPE-BASED METACOGNITION
% linear fit of all correct trials as function of Wagers of one subject
%         p_conf_post      = polyfit(1:n_wagers,PostCertain.PercSuc,1);
% 		r_conf_post      = p_conf_post(1) .* [1:n_wagers] + p_conf_post(2);
% 		p_err_post       = polyfit(1:n_wagers,PostCertain.PercUnsuc,1);
% 		r_err_post       = p_err_post(1) .* [1:n_wagers] + p_err_post(2);
% 		
%         p_conf_pre       = polyfit(1:n_wagers,PreCertain.PercSuc,1);
% 		r_conf_pre       = p_conf_pre(1) .* [1:n_wagers] + p_conf_pre(2);
% 		p_err_pre        = polyfit(1:n_wagers,PreCertain.PercUnsuc,1);
% 		r_err_pre        = p_err_pre(1) .* [1:n_wagers] + p_err_pre(2);
%         r_conf_corrected = r_conf_post(:) - r_conf_pre(:);
% 		r_err_corrected  = r_err_post(:)  - r_err_pre(:);
% 		
% 		meta_conf          = p_conf_post(1) - p_conf_pre(1);
% 		meta_err           = -p_err_post(1) + p_err_pre(1);
% 		meta_conf_pre      = p_conf_post(1) - p_conf_pre(1);
% 		meta_err_pre       = -p_err_post(1) + p_err_pre(1);
% 		General(file_index).meta_slope               = meta_conf(:) + meta_err(:);
% 		meta_pre           = meta_conf_pre + meta_err_pre;


Tab_Diff =  struct2table(TDiff);
Tab_Session_Diff = [Tab_Session_Diff;Tab_Diff];
end       


Tab_Session = struct2table(General);
save([path_save,filesep,  monkey, '\',monkey(1:3), '_VariablesPerSession_M2S_Wagering_' ,dateStart,'_until_', dateEnd] ,'Tab_Session');
writetable(Tab_Session,[path_save,filesep,  monkey, '\',monkey(1:3), '_VariablesPerSession_M2S_Wagering_' ,dateStart,'_until_', dateEnd,'.txt'], 'Delimiter', ',')
writetable(Tab_Session,['C:\Users\kkaduk\Dropbox\promotion\Projects\Wagering_Monkey\Data\'  monkey(1:3),filesep, monkey(1:3), '_VariablesPerSession_M2S_Wagering_' ,dateStart,'_until_', dateEnd,'.txt'], 'Delimiter', ',')

writetable(Tab_Session_Diff,[path_save,filesep,  monkey, '\',monkey(1:3), '_VariablesPerSession_ForEachDifficultyLevel_M2S_Wagering_' ,dateStart,'_until_', dateEnd,'.txt'], 'Delimiter', ',')
writetable(Tab_Session_Diff,['C:\Users\kkaduk\Dropbox\promotion\Projects\Wagering_Monkey\Data\'  monkey(1:3),filesep, monkey(1:3), '_VariablesPerSession_ForEachDifficultyLevel_M2S_Wagering_' ,dateStart,'_until_', dateEnd,'.txt'], 'Delimiter', ',')


disp(['saved  C:\Users\kkaduk\Dropbox\promotion\Projects\Wagering_Monkey\Data\' monkey(1:3),filesep, monkey(1:3),  '_VariablesPerSession_M2S_Wagering_' ,dateStart,'_until_', dateEnd,'.txt'])
