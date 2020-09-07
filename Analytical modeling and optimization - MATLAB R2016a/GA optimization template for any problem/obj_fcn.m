function [obj, c, ceq] = obj_fcn(x)

% this function designs an axial flux pcb motor

global results_array var_names counter 

%% constants

n_rpm = 2000;  %rpm, rated speeed
t_el = 0.3;    %Nm, nominal torque


%% optimization parameters

rout_mag = x(1)/2*1e-3;   %m, magnet outer radius
l_mag = x(2)/2*1e-3;  %m, magnet thickness
pole = x(3)*2;  % number of poles
slots = x(4);  %number of slots



%% basic relations



%% magnetic loading calculation







%% mass calculation


m_total = 1;  %kg, total mass of the motor

%% objective function

obj = m_total;  %kg, total mass of the motor



%% optimization constraints


% When there are integer constraints, ga does not accept linear or nonlinear equality constraints, only inequality constraints.
ceq = []; 

t_min = 0;
% inequality constraints defined below should be negative, buraya eklenen constraintler yukarýdaki exit if'lerine de eklenmeli NaN olarak
c(1) = t_min - t_el;  % air gap flux density should be less than bremenant constraint




%% print the results

if ( all ( c<=0 )  )  

results_array(counter,:) = [rout_mag*1e3 l_mag*1e3 pole slots m_total obj]; %kaydedilen sonuclar
var_names = {'rout_mag', 'l_mag' ,'pole', 'slots','m_total' , 'obj'};  % Variable names, get from objective function

counter = counter+1;

end


end