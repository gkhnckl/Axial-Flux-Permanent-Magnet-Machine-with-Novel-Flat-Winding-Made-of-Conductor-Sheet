function [ BIx_tot, BIy_tot, BIIx_tot, BIIy_tot,BIy_fund_peak] = FieldAnalysis_DSAFPM( hm, L, tau_p, br25, alpha_p, mur, mu0, N, RTC, temp_amb, x, y )

% No-load magnetic field analysis of a single side axial flux PM motor
% Source: "Field Computation for an Axial Flux Permanent-Magnet 
% Synchronous Generator" by T.F. Chan, 2009


n = 1:1:N; % Number of harmonics to be considered

%% Temperature dependency of Brem

br = br25 * ( 1 + RTC/100 * (temp_amb - 25) );  %T, Residual magnetism @temp_amb
% RTC has negative value and it is in %/C

%% Define magnetization coefficients for My

MnN = 4*br./(mu0*n*pi).*(sin(n*pi*alpha_p/2)).*sind(180*n/2).^2;  %Fourier series coefficients of magnetization vector in y-dir

%% Define variables in the matrix
E1 = diag(exp(n*pi*L/tau_p),0);
E2 = diag(exp(-n*pi*L/tau_p),0);
E3 = diag(exp(n*pi*hm/tau_p),0);
E4 = diag(exp(-n*pi*hm/tau_p),0);

%% Define E and C matrix

E = [E1 E2 zeros(N,N) zeros(N,N);  
zeros(N,N) zeros(N,N) diag(ones(1,N)) diag(ones(1,N));
E3 E4 -E3 -E4;
diag(-pi*n/tau_p,0)*E3 diag(pi*n/tau_p,0)*E4 diag(pi*n*mur/tau_p,0)*E3 diag(-pi*n*mur/tau_p,0)*E4;];

Y = [zeros(N,1);
zeros(N,1);
zeros(N,1);
MnN';];

%% Calculate coefficients

% C = inv(E)*Y;    %calculation of constants
C = E\Y;  %Calculation of constants
C1 = C(1:N); C2 = C(N+1:2*N);
C3 = C(2*N+1:3*N); C4 = C(3*N+1:4*N);

%% Calculate the field by lower and upper magnets for each harmonic

for i=1:N

%field by lower magnets
BIx_low(i,:)= pi*i*mu0/tau_p .* (C1(i).*exp(pi*i*y/tau_p)+C2(i).*exp(-pi*i*y/tau_p)) .* sin(pi*i.*x/tau_p);  %T, Bx in Reg1
BIy_low(i,:)= -pi*i*mu0/tau_p .* (C1(i).*exp(pi*i*y/tau_p)-C2(i).*exp(-pi*i*y/tau_p)) .* cos(pi*i.*x/tau_p);  %T, By in Reg1

BIIx_low(i,:)= pi*i*mu0*mur/tau_p .* (C3(i).*exp(pi*i*y/tau_p)+C4(i).*exp(-pi*i*y/tau_p)) .* sin(pi*i.*x/tau_p);  %T, Bx in Reg2
BIIy_low(i,:)= -pi*i*mu0*mur/tau_p .* (C3(i).*exp(pi*i*y/tau_p)-C4(i).*exp(-pi*i*y/tau_p)) .* cos(pi*i.*x/tau_p)...
    + mu0* MnN(i)*cos(pi*i*x/tau_p)  ;  %T, By in Reg2

%field by upper magnets
BIx_up(i,:)= pi*i*mu0/tau_p .* (C1(i).*exp(pi*i*(y+L)/tau_p)+C2(i).*exp(-pi*i*(y+L)/tau_p)) .* sin(pi*i.*x/tau_p);  %T, Bx in Reg1
BIy_up(i,:)= -pi*i*mu0/tau_p .* (C1(i).*exp(pi*i*(y+L)/tau_p)-C2(i).*exp(-pi*i*(y+L)/tau_p)) .* cos(pi*i.*x/tau_p);  %T, By in Reg1

BIIx_up(i,:)= pi*i*mu0*mur/tau_p .* (C3(i).*exp(pi*i*(y+L)/tau_p)+C4(i).*exp(-pi*i*(y+L)/tau_p)) .* sin(pi*i.*x/tau_p);  %T, Bx in Reg2
BIIy_up(i,:)= -pi*i*mu0*mur/tau_p .* (C3(i).*exp(pi*i*(y+L)/tau_p)-C4(i).*exp(-pi*i*(y+L)/tau_p)) .* cos(pi*i.*x/tau_p)...
    + mu0* MnN(i)*cos(pi*i*x/tau_p)  ;  %T, By in Reg2


end


%% Calculate field caused by both magnets using superposition

BIx = BIx_low + BIx_up;      %T, Bx in Reg1
BIy = BIy_low + BIy_up;     %T, By in Reg1
BIIx = BIIx_low + BIIx_up;      %T, Bx in Reg2
BIIy = BIIy_low + BIIy_up;       %T, By in Reg2

BIy_fund_peak = max(BIy(1,:)); %T, peak fundamental flux density at the air gap

%% Sum all harmonics

BIx_tot = sum(BIx);  %finding magnetic field by summing all harmonics
BIy_tot = sum(BIy);
BIIx_tot = sum(BIIx);
BIIy_tot = sum(BIIy); 


end

