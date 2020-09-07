function t_core = deflection(rin_mag, rout_mag, rin_fix, ya_max, alpha_p, b_peak)

% this file calculates cimcumferential deflection of a circular disk due to
% maxwell stress
% reference: roark's formulas for stress and srain
% pdf page numbers : 464 and 472



%% constants

a = rout_mag;     %m, magnet outer radius
b = rin_fix;        %m, shaft radius
r0 = rin_mag;         %m, magnet inner radius
% t = t_core;             %m, core thickness

E = 200e9;   %Pa, young's modulus of iron = 210e9
v = 0.29;   %no unit, poissons ratio of iron = 0.29
u0 = pi*4e-7;   %H/m, permeability of free space

% flat circular plate constants
c2 = 1/4 * ( 1-(b/a)^2*(1+2*log(a/b)) );
c3 = b/(4*a) * (  ((b/a)^2+1)*log(a/b) + (b/a)^2 -1  );
c8 = 1/2 * (  1 + v + (1-v)*(b/a)^2  );
c9 = b/a * (  (1+v)/2*log(a/b) + (1-v)/4*(1-(b/a)^2)  );

% flat circular plate loading constants
L11 = 1/64 * (   1 + 4*(r0/a)^2 - 5*(r0/a)^4 - 4*(r0/a)^2*( 2+(r0/a)^2 )*log(a/r0) );
L17 = 1/4 * (   1 - (1-v)/4*(1-(r0/a)^4) - (r0/a)^2*(1+(1+v)*log(a/r0))   );


%% deflection calculation

q_actual = b_peak^2/(2*u0);    %Pa or N/m2,   maxwell stress tensor
q = q_actual * alpha_p;    %Pa or N/m2, stress in tüm diske yayýlan ortalama degeri


% D = E*t^3 / (12*(1-v^2));    %Nm, plate constant
Qb = q/(2*b) * (a^2-r0^2);   %N/m, unit shear force
Mrb = -q*a^2/c8 * (  c9/(2*a*b)*(a^2-r0^2) - L17  );  %N, unit radial bending moment at the shaft

% ya = Mrb*a^2*c2/D + Qb*a^3*c3/D - q*a^4*L11/D;  %m, maximum deflection occuring at tip of core
% ya = abs(ya);

t_core = (   12*a^2*(1-v^2)/(E*ya_max) * (Mrb*c2 + Qb*a*c3 - q*a^2*L11)   )^(1/3);   %m, calculate reqired thickness from max ya
t_core = abs(t_core);

end