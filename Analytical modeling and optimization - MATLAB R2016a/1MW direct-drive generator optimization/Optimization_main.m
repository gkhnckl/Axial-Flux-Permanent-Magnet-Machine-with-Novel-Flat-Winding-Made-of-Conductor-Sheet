clear all; clc;

global counter results_array var_names
counter = 1;
results = [];


%% borders

pole_max = 100;   % number of poles
pole_min = 10;  

dia_out_max = 4;   %m, outer diameter  
dia_out_min = 0.1;   

width_fw_max = 30e-3;   %m, flat wire width
width_fw_min = 5e-3;   

thick_fw_max = 5e-3;  %m, flat wire thickness
thick_fw_min = 1e-3;

J_rms_max = 10e6;  %A/m2, rms current density
J_rms_min = 1e6;

stage_max = 8;  % number of stages for multi-stage optimization
stage_min = 1;

paral_max = 14;
paral_min = 1;

fw_material_max = 1;
fw_material_min = 1;

UB = [pole_max/2 dia_out_max*100 width_fw_max*1e4 thick_fw_max*1e3 J_rms_max*1e-5 stage_max paral_max fw_material_max];   % Boundaries for multi-stage design
LB = [pole_min/2 dia_out_min*100 width_fw_min*1e4 thick_fw_min*1e3 J_rms_min*1e-5 stage_min paral_min fw_material_min];


%% ga options 

tic
options = gaoptimset(@ga);
% options.TolFun = 0.1;    %obj fcn daki change ne kadar önemli, 1e-4
% options.StallGenLimit = 100;  %if the average change in the fitness function value over Stall generations is less than Function tolerance, the algorithm stops. 200
% options.CrossoverFraction = 0.95;  %Elit sayýsý haricinde kalanlarýn kaçý mutasyon ile next generasyona aktarýlacak, 0.8
% options.PopulationSize = 400;  % nufüs,Positive integer | {50} when numberOfVariables <= 5, {200} otherwise | {min(max(10*nvars,40),100)} for mixed-integer problems,400
% options.TolCon = 1e-1;  % constraint tolerance, 1e-4
% options.EliteCount = 5;  % number of elites


% options.Display ='iter';  % show the results on commend windows at each iteration 
options.PlotFcns = {@gaplotbestf};  % plot the required plots https://www.mathworks.com/help/gads/genetic-algorithm-options.html#f14474
options.PlotInterval = 1 ;  %defines plot period in terms of generation numbers
% Best fitness ('gaplotbestf') plots the best score value and mean score versus generation.
% Stopping ('gaplotstopping') plots stopping criteria levels.
% Distance ('gaplotdistance') plots the average distance between individuals at each generation.
% Range ('gaplotrange') plots the minimum, maximum, and mean score values in each generation.



% [x,fval,exitflag,output,population,scores] = ga(@obj_fcn,numel(LB),[],[],[],[],LB,UB,@nonlcon,[1,2,3,4,5],options);  % optimization for single-stage
[x,fval,exitflag,output,population,scores] = ga(@obj_fcn_MS,numel(LB),[],[],[],[],LB,UB,@nonlcon,[1,2,3,4,5,6,7,8],options);  % optimization for multi stage design


results_array = unique(results_array,'rows');
results_array = sortrows(results_array,size (results_array,2)); %obj fcn a göre sýrala

results_Table = array2table(results_array,'VariableNames',var_names);
toc