%Function to calculate structural mass of electrical machines

%   INPUTS
% Torque_per_stage: Mechanical torque per stage, Nm
% r_out_mag: Outer radius of the magnets, m
% r_in_mag:  Inner radius of the magnets, m
% endwind_in_fw: End winding length of the flat wire at inner radius, m
% b_peak: Peak flux density at the air gap, T
% embrace: Emrace of the magnets
% g_cl:  Air gap clearance, m
% dist_rotor: Distance between two rotor disks, m
% t_core_sat: Rotor back core thickness calculated by saturation
% stage: number of stages

%    OUTPUTS
% mass_rotor_structure: Total structural mass of rotor (!Mass of shaft should be added)
% mass_stator_structure: Total structural mass of stator (!Mass of bearing may be added)

function [mass_rotor_structure, mass_stator_structure, m_rotor_core, t_core_in, t_core_out] = structure_MS(Torque_per_stage, r_out_mag, r_in_mag, endwind_in_fw, b_peak,embrace, g_cl, dist_rotor,t_core_sat, stage)


n = 5;          % Number of Structural Arms for the Rotor
rad_cl = 15e-3;  %m, radial clearance between rotor and stator

%% Constants %%

E = 2e11;          % Young's Modulus (Pa)
rho_steel = 7800;  % Density of Steel (kg/m^3)


r_shaft = r_out_mag/8;         % Shaft Radius
hr = max(15e-3, r_out_mag / 100);  %m, rotor structural yoke thickness
hs = max(15e-3, r_out_mag / 100);  %m, stator structural yoke thickness
r_in_stator = r_in_mag - endwind_in_fw;  %m, stator inner radius



%% axial deflection due to maxwell stress to determine core thickness
% this deflection is due to the attraction force between two rotors

percent_def_max_ax = 8;  % percent, allowed deflection in percentage of airgap clearance
def_max_ax = g_cl * percent_def_max_ax/100;  %m, allowed deflection in one side of the rotor

t_core_def = deflection(r_in_mag, r_out_mag, (r_in_stator-rad_cl), def_max_ax, embrace, b_peak); %m, minimum core thickness according to deflection criteria



%% circumferential deflection due to shear stress
% this deflection is due to mechanical torque of the motor and determines
% rotor arm shape of the structure

sigma = 2*Torque_per_stage*stage / ( pi*(r_out_mag-r_in_mag)*(r_out_mag+r_in_mag)^2 );  %Pa, N/m2, Sheer Stress (Pa)

% b = 2*r_out_mag*r_shaft / n;        %m, hallow rotor arm genisligi
b=2*r_shaft*sin(2*pi/n/2); %m, hallow rotor arm genisligi
d = r_out_mag / 100;          %m, hallow rotor arm axial length
t_w = r_out_mag / 100;        % hallow rotor arm et kalýnlýðý

I_arm_tor = ((d.*b.^3)-((d-2.*t_w).*(b-2.*t_w).^3))./12;         %m^4, Second Moment of the Area of the Rotor's Arm
l_i = ( r_in_stator - rad_cl - hr ) -  r_shaft;                  %m, Length of Rotor Arm (m)

z = pi* (r_out_mag^2-r_in_mag^2) * sigma.*(l_i.^3)./n/3/E/I_arm_tor;             %m, Circumferential Deflection of the Rotor Arm (m)

z_all = 0.005.*2.*pi.*(r_out_mag)./360;       % Allowable Circumferential Deflection of the Rotor Arm (m)

while z_all < z
        d = d+0.002;
        I_arm_tor = ((d.*b.^3)-((d-2.*t_w).*(b-2.*t_w).^3))./12;
        l_i = ( r_in_stator - rad_cl - hr ) -  r_shaft;                  %m, Length of Rotor Arm (m)
        z = pi* (r_out_mag^2-r_in_mag^2) * sigma.*(l_i.^3)./n/3/E/I_arm_tor;             %m, Circumferential Deflection of the Rotor Arm (m)
end
 
a = (b.*d)-((b-2.*t_w).*(d-2.*t_w));  %m2, cross sectional area of the rotor arm


%% Mass calculation

t_core_out = max(t_core_sat, t_core_def);  %m, outer core thickness
t_core_in = max(15e-3,r_out_mag/150);    %m, inner core thickness
len_st = stage * dist_rotor + (stage-1)*t_core_in;  %m, axial length of the stator structure


%rotor structural mass
r_out_rotor = r_out_mag;  %m, outer radius of the rotor core
r_mid_rotor = r_in_stator - rad_cl;  %m, middle radius of the rotor
r_in_rotor = r_mid_rotor - hr; %m, inner radius of the rotor

m_rotor_core1 = pi*(r_out_rotor^2-r_mid_rotor^2)* (2*t_core_out + (stage-1)*t_core_in) * rho_steel;  %kg, core mass of double sided rotor back core
m_rotor_core2 = pi*(r_mid_rotor^2-r_in_rotor^2) * (len_st+2*t_core_out) * rho_steel; %kg, core mass for rotor plate

m_rotor_core = m_rotor_core1 + m_rotor_core2; %kg, rotor core mass, note that this belongs to the active mass,not structurral mass

mass_rotor_structure = (n.*l_i.*a.*rho_steel);        %kg, Rotor mass (includes only rotor arm)

%stator structural mass
r_out_s = r_out_mag + rad_cl + hs; %m, outer radius of stator structure
r_in_s = r_out_mag + rad_cl ;  %m, inner radius of stator structure
mass_stator_structure = pi * (r_out_s^2 - r_in_s^2) * len_st * rho_steel;  %kg, stator structural mass




%note that rotor structural mass includes only rotor hallow arms. rotor
%core itself is included in the active mass


%% Print Results %%
%fprintf('MASS:\nMass of rotor structure = %d kg\nMass of stator structure = %d kg\nTotal Mass = %d kg\n\n',round(mass_rotor_structure), round(mass_stator_structure), round(total_mass))

end



