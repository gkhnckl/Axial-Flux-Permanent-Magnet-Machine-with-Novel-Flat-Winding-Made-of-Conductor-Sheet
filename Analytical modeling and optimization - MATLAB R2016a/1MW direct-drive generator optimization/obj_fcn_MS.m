function [fitness, c, ceq] = obj_fcn_MS(x)
% clc; clear all;
% Fitness function for multi-stage design
global results_array counter var_names

%% optimization parameters


pole = x(1)*2;   % number of poles
dia_out = x(2)/100;  %m, outer radius
width_fw = x(3)/1e4;  %m, flat wire width
thick_fw = x(4)/1e3;  %m, flat wire thickness
J_rms = x(5)/1e-5;  %A/m2, rms current density
stage = x(6); % number of stages
paral = x(7);   % number of parallel turns per phase (pole pairs can be connected in parallel or series)
fw_material = x(8); % flat wire material


% pole = 46;
% dia_out = 4.77;
% width_fw = 18.6e-3;
% thick_fw = 3e-3;
% J_rms = 3.8e6;
% stage = 2;
% paral = 2;

%% constants

p_out = 5e6;  %W, rated output power
n_rpm = 400;  %rpm, rated speed of the machine

p_out_stage = p_out / stage;  %W, output power
m = 3;  % phase number
thick_iso = 0.18e-3;  %m, isolation paper thickness
dist_fw_iso = 0.02e-3; %m, clearance between isolation paper and flat wire
dist_fw_mid = 0.5e-3; %m, flat wire middle space distance, U'nýn arasý
endwind_in_fw = 2*width_fw;  %m, flat wire end winding length at inner radius
endwind_out_fw = 15e-3;  %m, flat wire end winding length at outer side, welding place

temp_amb = 80;  %deg, operating temperature
ur = 1.05;  %relative permeabity of magnets
b_sat = 1.7;  %T, core saturation flux density
mu0 = pi*4e-7; %air permeability
alpha_p = 0.78; %embrace of magnets 0.78 approximation

g_cl = round(dia_out * 10)/1e4;   %m, air gap clearence between flat wire and magnets, it is 0.1% of outer diameter, ref: mueller'in makalesi


% flat wire material properties
resistivity = [1.68e-8 2.65e-8];  %ohm*m, resistivity of flat wire material at 20 deg
temp_coef = [0.00386 0.004229];  % temperature dependency of resistivity
density_cond = [8.96e3 2.7e3]; %kg/m3, density of flat wire material conductor
sp_cost_cond = [15 4];  %€/kg, specific cost of the flat wire material




%% frequency calculation 

frequ = pole*n_rpm/120;  %Hz, back emf frequency
t_per = 1/frequ;  %sec, period


%% maximum flat wire calculation

r_out = dia_out / 2;  %m, outer radius of the core
r_in = r_out * (1-tan(pi/pole))/(1+tan(pi/pole));  %m, inner radius of the core
r_mean = (r_out+r_in)/2;  %m, mean radius

% limitation at the inner radius, num_fw_max_in
circ_fw_in = 2*pi*(r_in-endwind_in_fw);  %m, available circumferential distance at the inner radius
num_fw_max_in = floor (circ_fw_in / ( thick_fw+thick_iso+dist_fw_iso ));  %max number of fw that can be used according to inner side limit 

% limitation at the bending point at the middle, num_fw_max_mid
r_fw_bend = ( r_out-r_in ) / 2 / sind(180/pole) + thick_fw ;  %m, radius at which fw bend is occured at the middle
circ_fw_mid = 2*pi * r_fw_bend;  %m, available circumferential distance at the bend radius
num_fw_max_mid = floor( circ_fw_mid / ( (thick_fw+thick_iso)/sind(45-180/pole) ) );  %max mumber of fw that can be used according to bending limit

num_fw_max = min( num_fw_max_in , num_fw_max_mid);  %max number of fw allowed

%% total number of flat wires calculation

% num_fw, 3'ün katý olmalý
% iki fw arasýndaki açý 60'ýn böleni olmalý

num_fw_init = num_fw_max;  % max number of fw allowed geometrically
alpha = 180*pole / num_fw_max;  %deg, electrical angle between two fws


if( alpha > 60 ) % if initial alpha higher than 60, skip the iteration
    k_alpha_pen = 1e8;  % penalty multiplier for angle between two fws
    fitness =  (alpha-60) * k_alpha_pen;  % if design is ridiculus, give a large penalty
    c(1) = alpha-60;  % say that constraints are also not satisfied
    c(2)=NaN; c(3)=NaN; c(4)=NaN;  % say other constraints are also not satisfied
    ceq = [];
    return
end
   
check1 = logical( rem(num_fw_max,m) );  % num_fw 3'ün katý mý?
check2 = logical ( rem( 60,alpha ) ); % iki fw arasýndaki açý 60'ýn böleni mi?

while (check1 || check2) % if one of the criterias are not satisfied, decrease the number of total fws
   num_fw_max = num_fw_max - 1;  %kriterler saðlanmýyorsa bir eksilt   
   alpha = 180*pole / num_fw_max;  %deg, new electrical angle between two fws
   
   check1 = logical( rem(num_fw_max,3) );  % num_fw 3'ün katý mý?
   check2 = logical ( rem( 60,alpha ) ); % iki fw arasýndaki açý 60'ýn böleni mi?
 
end % eksiltmeye iki kriter saglanana kadar devam et

num_fw = num_fw_max;  % total number of flat wires in stator per stage

stator_ff = num_fw / num_fw_init * 100;  % how much allowable space used?



%% Phase current calculation

area_fw = thick_fw * width_fw;  %m2;  cross sectional area of a fw.
I_ph_stage_rms = area_fw * J_rms * paral;  %A;  phase current per stage considering number of parallel turns per phase



%% Phase voltage calculation

V_ph_fund_stage_rms = p_out_stage / ( m * I_ph_stage_rms);  %V, phase voltage per stage


%% winding distribution factor calculation, it is defined over one pole pair

num_fw_phase = num_fw / m;  % total number of flat wires used per phase per stage
num_fw_phase_series = num_fw_phase / paral;   % number of flat wires in series used per phase per stage
num_fw_loop = pole / 2; % number of flat wires in a loop
num_loop_ph = num_fw_phase / num_fw_loop;  % number of loops in a phase assuming all connected in series

alpha = 180 * pole / num_fw;  %deg, electrical angle between two zigzag loops
q = 60 / alpha;  % number of loops per phase segment (A -C B -A C -B)
k_d = sind(q*alpha/2) / (q * sind(alpha/2)); %distribution factor



%% B_rms_fund calculation 

edge_mag = (r_out-r_in)/sqrt(2);  %m, square magnet edge length
len_cond = edge_mag;  %m, conductor length in a fw. there are four conds in a fw

r_fw_mid_bot = sqrt( ( len_cond /(2*sqrt(2)) )^2 + ( r_in+len_cond/(2*sqrt(2)) )^2 );  %m, mean radius of bottom conductors in a fw
r_fw_mid_top = sqrt( ( len_cond /(2*sqrt(2)) )^2 + ( r_out-len_cond/(2*sqrt(2)) )^2 );  %m, mean radius of top conductors in a fw
r_fw_mid = mean([r_fw_mid_bot r_fw_mid_top]);  %m, mean radius of all conductors in a fw.


b_rms_fund = V_ph_fund_stage_rms * 30*sqrt(2) / ( pi * n_rpm * r_fw_mid * len_cond * 2*num_fw_phase_series * k_d );  %T, rms flux density of fundamental component



%% Temperature dependency of Brem
% Assuming N42M magnet is used. It has max temp of 100C

br25 = 1.29; %T, residual magnetism at 25C
RTC = -0.11;  % Percent/C, reversible temperature coefficient of residual induction
br = br25 * (1+(RTC/100)*(temp_amb-25));  %T, residual magnetism at ambient temp


%% Magnet thickness calculation

b_peak_fund = b_rms_fund * sqrt(2);   %T,peak fundamental flux density at the air gap
b_avg_fund =  b_peak_fund * 2 / pi;  %T, average fundamental flux density at the air gap


if( b_peak_fund > br ) % quit for peak flux density higher than remenance flux density
 
    k_br_pen = 1e9; % penalty multiplier for flux density higher than Br
    fitness = k_br_pen * (b_peak_fund - br);
    c(1) = b_peak_fund - br; % say that constraints are also not satisfied
    c(2)=NaN; c(3)=NaN; c(4)=NaN;  % say other constraints are also not satisfied
    ceq = [];
    return 
end

% % initial value of magnet thickness with primitive methods
g_magnetic = 2*g_cl+2*width_fw+dist_fw_mid;  %m, magnetic air gap from magnet to magnet
t_mag_init = (b_peak_fund) * ur * g_magnetic / (2*(br-b_peak_fund));    %m, magnet thickness


% % Magnet thickness calculation
harm = 15;  % harmonics to be considered
g_field = g_magnetic + 2*t_mag_init;  %m, distance between core faces
tau_p  = (2*pi*r_mean) / pole;  %m, pole pitch in circumferential direction
x_field = linspace(0.5*tau_p, 2.5*tau_p, 1001);    %m, circumferential length on which field is investigated
y_field = g_field/2;   %m, axial distance on which field is investigated

b_peak_fund_aim = b_peak_fund;   %T, aimed fundamental air gap peak flux density
t_mag_cur = t_mag_init;  %m, current value of magnet thickness
[~,~,~,~,b_peak_fund_cur] = FieldAnalysis_DSAFPM( t_mag_cur, g_field, tau_p, br25, alpha_p, ur, mu0, harm, RTC, temp_amb, x_field, y_field );

error_B = abs(b_peak_fund_cur-b_peak_fund_aim);  %T, error in fundamental air gap peak flux density
error_ref_B = 0.005;  %T, allowed maximum error in B


while (error_B > error_ref_B)

    t_mag_cur = t_mag_cur * ( 2 - b_peak_fund_cur/b_peak_fund_aim); %m, updated magnet thickness
    g_field = g_magnetic + 2*t_mag_cur;  %m, distance between core faces
    y_field = g_field/2;   %m, axial distance on which field is investigated
    [~,~,~,~,b_peak_fund_cur] = FieldAnalysis_DSAFPM( t_mag_cur, g_field, tau_p, br25, alpha_p, ur, mu0, harm, RTC, temp_amb, x_field, y_field );
    error_B = abs(b_peak_fund_cur-b_peak_fund_aim);  %T, error in fundamental air gap peak flux density
   
end


% % Update final value with rounded magnet thickness

t_mag = round(t_mag_cur,4);  %m, final value of magnet thickness

g_field = g_magnetic + 2*t_mag;  %m, distance between core faces
y_field = g_field/2;   %m, axial distance on which field is investigated
[~,BIy_tot,~,~,b_peak_fund] = FieldAnalysis_DSAFPM( t_mag, g_field, tau_p, br25, alpha_p, ur, mu0, harm, RTC, temp_amb, x_field, y_field );


%% Torque calculation

Torque_per_stage = m * b_rms_fund * r_fw_mid * len_cond/sqrt(2) * 2*num_fw_phase_series * k_d * I_ph_stage_rms;  %Nm, rated torque of the generator


%% phase resistance calculation of active flat wires

ro_cu_20 = resistivity(fw_material);  %ohm*m, resistivity of the conductor at 20deg
c_temp = temp_coef(fw_material);  %material resistivity temperature coefficient
ro_cu = ro_cu_20 * (1 + c_temp * (temp_amb-20));  %ohm*m, resistivity of the conductor at operating temp


len_fw = 4*len_cond + 4*width_fw;   %m, conductor length of a single flat wire
len_ph = len_fw * num_fw_phase;  %m, conductor length for one phase assuming all fws are connected in series

res_ph_all_series = ro_cu * len_ph / area_fw;  %ohm, phase resistance per stage per phase assuming all pole pairs are connected in series (paral=1)
res_ph_act = res_ph_all_series / paral^2;  %ohm, phase resistance per phase per stage


%% phase resistance calculation of end windings

r_endwind = r_out + 2*endwind_out_fw;  %m, approximated radius of end windings
circ_endwind = 2 * pi * r_endwind;  %m, circumferential distance of one circle of end winding

if rem(paral,2) == 0 % if number of parallel turns per phase is even   
    len_endwind = (paral/2-1) / (paral/2) * 2 * circ_endwind * num_loop_ph;   %m, apprx. end winding length 
else % if number of parallel turns per phase is odd
    len_endwind = (paral-1) / paral * 2 * circ_endwind * num_loop_ph;   %m, apprx. end winding length
end

area_endw_mult = 2;  % cross sectional area of end winding divided by cross sectional area of fw
res_endw = ro_cu * len_endwind / (area_fw*area_endw_mult) / paral^2;  %ohms, resistance of end windings per stage

res_ph = res_ph_act + res_endw;  %ohm, phase resistance per stage


%% I2R losses

p_cu_loss = m * I_ph_stage_rms^2 * res_ph * stage;   %W, joule copper losses for all stages


%% eddy current loss calculation
count = 1;

for eddy_lamination = linspace(0,1,10);

    % in this calculation, a side of flat wire is divided into multiple
    % regions in axial direction and field is calculated at each regions.
    % The eddy currents are calculated using the field derivation at each
    % region and to find actual eddy currents, found eddy current losses
    % according to regional field calculations are averaged
    
%first calculate the field on conductors
y_field = g_field/2 + dist_fw_mid/2 + width_fw*eddy_lamination ;   %m, axial distance at the middle of conductor
[B_eddy_x,B_eddy_y] = FieldAnalysis_DSAFPM( t_mag, g_field, tau_p, br25, alpha_p, ur, mu0, harm, RTC, temp_amb, x_field, y_field );

% harmonics analysis for b field in circumferential (x) direction
h_max = 15;   % max number of harmonics to be considered
B_eddy_x_fft = resample(B_eddy_x,2048,numel(B_eddy_x)); % adjust array size to 2048 points
L = 2048;                 % fft calculation
NFFT = 2^nextpow2(L);     % fft calculation
fft_res = fft(B_eddy_x_fft, NFFT)/L; % fft calculation
B_eddy_x_harm_peak = 2*abs(fft_res(2:h_max+1))';   % peak value of each harmonic
% bar(B_eddy_x_harm_peak)   % harmonics are calculated in this step, bar gives peak value of each harmonic

% harmonics analysis for b field in axial (y) direction
% h_max = 15;   % max number of harmonics to be considered
B_eddy_y_fft = resample(B_eddy_y,2048,numel(B_eddy_y)); % adjust array size to 2048 points
L = 2048;                 % fft calculation
NFFT = 2^nextpow2(L);     % fft calculation
fft_res = fft(B_eddy_y_fft, NFFT)/L; % fft calculation
B_eddy_y_harm_peak = 2*abs(fft_res(2:h_max+1))';   % peak value of each harmonic
% bar(B_eddy_y_harm_peak)   % harmonics are calculated in this step, bar gives peak value of each harmonic

% ksp_x calculation, specific loss calculation corrector constant due to harmonics in the circumferential field
harm_eddy = 1:2:h_max ;  % odd harmonics considered in harmonic effect
ksp_x = sum(harm_eddy.*harm_eddy.*B_eddy_x_harm_peak(harm_eddy)'.*B_eddy_x_harm_peak(harm_eddy)') / max(B_eddy_x)^2;  % correction coefficient due to harmonics contribution to eddy losses

% ksp_y calculation, specific loss calculation corrector constant due to axial harmonics in the field
% harm_eddy = 1:2:h_max ;  % odd harmonics considered in harmonic effect
ksp_y = sum(harm_eddy.*harm_eddy.*B_eddy_y_harm_peak(harm_eddy)'.*B_eddy_y_harm_peak(harm_eddy)') / max(B_eddy_y)^2;  % correction coefficient due to harmonics contribution to eddy losses

% ksp2 calculation, specific loss calculation corrector due to tangential
% field at the middle of conductor

B_eddy_1 = B_eddy_y_harm_peak;  %T, magnetic field for eddy currents at axial direction
B_eddy_21 = B_eddy_x_harm_peak * cosd(45);  %T, magnetic field for eddy currents for bottom conductors at circumferential direction
B_eddy_22 = B_eddy_x_harm_peak * cosd(45+360/pole);  %T, magnetic field for eddy currents for top conductors at circumferential direction
B_eddy_31 = B_eddy_x_harm_peak * sind(45);  %T, magnetic field for eddy currents for bottom conductors at radial direction
B_eddy_32 = B_eddy_x_harm_peak * sind(45+360/pole);  %T, magnetic field for eddy currents for top conductors at radial direction

p_sp_1 = (thick_fw)^2 * (2*pi*frequ)^2 * B_eddy_1(1)^2 / 24 / ro_cu * ksp_y; % W/m3,  specific loss density of eddy current on the coils at axial direction
p_sp_21 = (width_fw)^2 * (2*pi*frequ)^2 * B_eddy_21(1)^2 / 24 / ro_cu * ksp_x; % W/m3,  specific loss density of eddy current on the bottom coils at circumferential direction
p_sp_22 = (width_fw)^2 * (2*pi*frequ)^2 * B_eddy_22(1)^2 / 24 / ro_cu * ksp_x; % W/m3,  specific loss density of eddy current on the top coils at circumferential direction
p_sp_2 = mean([p_sp_21 p_sp_22]);  %W/m3, specific loss density of eddy current on the coils at circumferential direction
p_sp_31 = (thick_fw)^2 * (2*pi*frequ)^2 * B_eddy_31(1)^2 / 24 / ro_cu * ksp_x; % W/m3,  specific loss density of eddy current on the bottom coils at radial direction
p_sp_32 = (thick_fw)^2 * (2*pi*frequ)^2 * B_eddy_32(1)^2 / 24 / ro_cu * ksp_x; % W/m3,  specific loss density of eddy current on the top coils at radial direction
p_sp_3 = mean([p_sp_31 p_sp_32]);  %W/m3, specific loss density of eddy current on the coils at radial direction
p_sp = p_sp_1 + p_sp_2 + p_sp_3;    %W/m3, total specific loss density of eddy current on the coils

vol_cu_eddy = thick_fw * width_fw * len_cond * 4 * num_fw; %m3, volume of the flat wires per stage
p_eddy_lamination(count) = p_sp * vol_cu_eddy * stage;  %W, eddy current loss on the coils for all stages per lamination of flat wire
count = count+1;
end

p_eddy = mean(p_eddy_lamination);  %W, eddy current loss on the coils for all stages

%% frictional loss and total power loss calculation

p_rot_bearing = p_out_stage * stage * (0.5/100);  %W, rotational and bearing losses are assumed to be 0.5% of the rated power
p_in = p_out_stage*stage + p_cu_loss + p_eddy + p_rot_bearing; %W, input power for all stages

eff = p_out_stage*stage / p_in * 100;  % efficiency in percent for all stages


%% outer core thickness according to saturation

flux_pp_peak = (r_out-r_in) * trapz(x_field , abs(BIy_tot)) / 2  ;   %Wb, peak flux per pole
t_core_sat = flux_pp_peak / (2*b_sat*(r_out-r_in));   %m, core thickness


%% structural mass calculation

dist_rotor = g_magnetic + 2*t_mag;  %m, distance between two rotor back cores
[mass_rotor_structure, mass_stator_structure, m_rotor_core, t_core_in, t_core_out] = structure_MS(Torque_per_stage, r_out, r_in, endwind_in_fw, max(BIy_tot),alpha_p, g_cl,dist_rotor,t_core_sat,stage);

m_structural = mass_rotor_structure + mass_stator_structure;  %kg, total structural mass for all stages

%% axial length calculation

len_axial = stage * dist_rotor + (stage-1) * t_core_in + 2*t_core_out;  %m, axial length of the generator

%% active and total mass calculation

dens_cu = density_cond(fw_material);     %kg/m3, density of the conductor
dens_mag = 7.5e3;     %kg/m3, density of magnets
dens_ins = 1.54e3;    %kg/m3;  density of expoxy insulation

vol_cu = (4*len_cond+endwind_in_fw+2*endwind_out_fw) * thick_fw *  width_fw * num_fw * stage; %m3, volume of the flat wires for all stages

m_mag = edge_mag^2*t_mag*pole*2*dens_mag*stage;   %kg, total magnet mass for all stages
m_cu = vol_cu * dens_cu;  %kg, total copper weight for all stages
m_ins = pi*(r_out^2-r_in^2)*thick_iso*dens_ins*stage;  %kg, insulation expoxy mass for all stages

m_active = m_rotor_core + m_mag + m_cu + m_ins;  %kg, total active mass


m_total = m_active + m_structural;  %kg, total mass including active and structural mass


%% cost estimation

sp_cost_cu = sp_cost_cond(fw_material);  % €/kg, specific cost of the flat wire material
sp_cost_mag = 60;  % €/kg, specific cost of the magnets
sp_cost_core = 2;  % €/kg, specific cost of the back iron
sp_cost_steel = 2;  % €/kg, specific cost of the structural steel



cost_mag = sp_cost_mag * m_mag;  % €, cost of the permanent magnets for all stages
cost_cu = sp_cost_cu * m_cu;  % €, cost of the copper for all stages
cost_core = sp_cost_core * m_rotor_core;  % €, cost of the core material for all stages
cost_steel = sp_cost_steel * m_structural; % €, cost of the structural steel for all stages

cost = cost_mag + cost_cu + cost_core + cost_steel; % €, total cost of the generator including active and inactive materials for all stages

%% fitness function


fitness = m_total;  % objective function
% fitness = cost;  % objective function


%% constraints

eff_min = 91;  % minimum efficiency in percent
V_ph_fund_stage_rms_max = 690/sqrt(3);  %V, Line to line voltage should be max 690V
% V_ph_fund_stage_rms_max = Inf;  %V, Line to line voltage should be max 690V

% When there are integer constraints, ga does not accept linear or nonlinear equality constraints, only inequality constraints.
ceq = []; 

% inequality constraints defined below should be negative, buraya eklenen constraintler yukarýdaki exit if'lerine de eklenmeli NaN olarak
c(1) = b_peak_fund - br;  % air gap flux density should be less than bremenant constraint
c(2) = eff_min - eff;   % minimum efficiency constraint
c(3) = V_ph_fund_stage_rms - V_ph_fund_stage_rms_max;   %V, voltage constraint, let it be 690 Vll
c(4) = abs( rem(pole,paral) );  % pole number should be multiplier of number of parallel turns per phase (paral)


%% if constraint are satisfied, record results

if ( all ( c<=0 )  )  

results_array(counter,:) = [pole dia_out width_fw*1e3 thick_fw*1e3 J_rms/1e6 stage paral fw_material frequ p_cu_loss/1e3 p_eddy/1e3 eff b_peak_fund stator_ff num_fw t_mag*1e3 edge_mag*1e3 len_axial Torque_per_stage*stage/1e3 I_ph_stage_rms V_ph_fund_stage_rms cost/1e3 m_mag m_active/1e3, m_structural/1e3, m_total/1e3, fitness]; %kaydedilen sonuclar
var_names = {'pole', 'Dia' ,'width_fw', 'thick_fw','J_rms' , 'Stage','Paral','FwMater','Freq', 'CopperLoss_kW', 'EddyLoss_kW', 'eff', 'b_peak_fund', 'StatorFF', 'NumFwSt', 'MagThick','MagEdge','AxialLength','Torque_kNm','I_ph_stage_rms', 'V_ph_fund_rms', 'Cost_kEuros', 'Mag_mass','ActMass', 'StrMass','TotalMass', 'Fitn'};  % Variable names, get from objective function

counter = counter + 1;

end

end