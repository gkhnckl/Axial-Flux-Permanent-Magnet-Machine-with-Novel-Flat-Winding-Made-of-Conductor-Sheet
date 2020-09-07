clear all
clc

% this function designs an axial flux pcb motor

global  results_array var_names counter

counter = 1;


%% borders

rout_mag_max = 70;   %mm, max outer radius
rout_mag_min = 35;   %mm, min outer radius

l_mag_max = 12;     %mm,  max magnet thickness
l_mag_min = 1;      %mm,  min magnet thickness

pole_max = 20;   %Arms, max number of poles
pole_min = 4;   %Arms, min number of poles

slots_max = 150;    %max number of slots
slots_min = 10;     %min number of slots

LB = [rout_mag_min*2 l_mag_min*2 pole_min/2 slots_min ];    % Boundaries for the optimization
UB = [rout_mag_max*2 l_mag_max*2 pole_max/2 slots_max ];


%% ga options 


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


tic
[x,fval,exitflag,output,population,scores] = ga(@obj_fcn,numel(LB),[],[],[],[],LB,UB,@nonlcon,[1,2,3,4],options);


results_array = unique(results_array,'rows');
results_array = sortrows(results_array,size (results_array,2)); %obj fcn a göre sýrala

results_Table = array2table(results_array,'VariableNames',var_names);
toc








