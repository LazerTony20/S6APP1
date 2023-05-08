% S6 APP1 E2023 - Calculs Problématique
% ROYA2019 - BLAD0901 - DESJ3602 
close all
clear
clc

% Constantes
q = 1.60*10^-19;    % (En coulomb)
e0 = 8.85*10^-14;   % Permittivité du vide (En F/cm)
Vt = 0.0259;        % En Volts (25.9mV)

% Valeurs connues
Vpp = 10;       % En Volts
Vmm = -10;      % En Volts
R1 = 50000;     % En Ohms 
R2 = 80000;     % En Ohms 
R3 = 10000000;  % En Ohms 
R4 = 120000;    % En Ohms 
R5 = 1;         % En Ohms 
R6 = 50;        % En Ohms
R7 = 220000;    % En Ohms 
R8 = 15000;     % En Ohms 
R10 = 0.1;      % En Ohms 
R11 = 200;      % En Ohms 

% Etage 1
Vt_E1 = Vt;
Beta_E1 = 200;  
Va_E1 = 100;        % En Volts
Idark = 5e-9;       % En Ampères
Vbe_approx = 0.7;   % En Volts
Rd = 1000000000;    % En Ohms

% Etage 2
k_E2 = 0.00303;
Vt_E2 = 0.7;
Va_E2 = 130;

% Etage 3
betaQ2      = 200;
betaQ3on    = 10;
betaQ3off   = 70;
betaQ3      = betaQ3off;
Vbe3_E3 = 0.64;
Vbe2_E3 = 0.59;
Vdel = 2.6;

%% Calculs ETAGE 1
disp('=====ETAGE1======')
Ve_E1 = 0 - Vbe_approx;
Vr3 = Ve_E1 - Vmm;
% Vr1 = Vpp - (Vmm + Vr3 + Vce_approx);
% Ic = Vr1/R1;
R3llRd = 1/((1/R3)+(1/Rd));
Ir3rd = Vr3/R3llRd;
Ie_E1 = Ir3rd + Idark;
Ic_E1 = Ie_E1*(Beta_E1/(Beta_E1+1));
Vc_E1 = Vpp - R1*Ic_E1;
Vceq_E1 = Vc_E1 - Ve_E1;
Ib_E1 = Ic_E1/(Beta_E1);
gm_E1 = Ic_E1/Vt_E1;
Rpi_E1 = Beta_E1/gm_E1;
R0_E1 = Va_E1/(Ic_E1);
disp(['Le courant Icq = ', num2str(Ic_E1*1000000000), ' nA'])
disp(['La tension Vceq = ', num2str(Vceq_E1), 'V'])
disp(['gm = ', num2str(gm_E1*1000), ' mSiemens'])
disp(['Rπ = ', num2str(Rpi_E1/1000000), ' MΩ'])
disp(['R0 = ', num2str(R0_E1/1000000), ' MΩ'])
disp(' ')
%% Calculs ETAGE 2
disp('=====ETAGE2======')
% Vov = [Vovx - Id]
Vovx_E2 = Vt_E2 - 8;
Ids = roots([1 (2/k_E2 - 2*Vovx_E2) (Vovx_E2).^2]);
Id_E2 = Ids(2);
Vsg_E2 = -(8 + Id_E2);
Vov_E2 = Vt_E2 - abs(Vsg_E2);
Vsd_E2 = -(Vpp-Vmm) - Id_E2*(R6+R5);
gm_E2 = (2*Id_E2)/Vov_E2;
R0_E2 = Va_E2/abs(Id_E2);   % P.390
disp(['Le courant Idq = ', num2str(abs(Id_E2)*1000), ' mA'])
disp(['La tension Vovq = ', num2str(Vov_E2), ' V'])
disp(['La tension Vdsq = ', num2str(Vsd_E2), ' V'])
disp(['gm = ', num2str(gm_E2*1000), ' mSiemens'])
disp(['R0 = ', num2str(R0_E2/1000), ' kΩ'])
disp(' ')
%% Calculs ETAGE 3
disp('=====ETAGE3======')
Vth_E3 = Vmm + (Vpp-Vmm)*R8/(R8+R7);
Rth_E3 = (R8.*R7)./(R8+R7);

Ib3_E3 = (Vth_E3 - Vmm - Vbe3_E3 - Vbe2_E3)/Rth_E3;
Ic2_E3 = betaQ2*(betaQ3+1)*Ib3_E3;
Ie2_E3 = Ic2_E3*(betaQ3+1)/betaQ3;
Vce2_E3 = Vpp - Vmm - Vdel - Ie2_E3*R10;
Vce3_E3 = Vpp - Vbe2_E3 - Vmm - Vdel - Ie2_E3*R10;

Ib2_E3 = Ie2_E3 - Ic2_E3;
Ie3_E3 = Ib2_E3;
Ic3_E3 = Ie3_E3*(betaQ3/(betaQ3 + 1));

disp(['La tension de l`équivalent de Thévenin à la base de Q3 = ', num2str(Vth_E3), ' V'])
disp(['La résistance de l`équivalent de Thévenin à la base de Q3 = ', num2str(Rth_E3/1000), ' kΩ'])
disp(['Le courant IBQ de Q3 = ', num2str(Ib3_E3*1000000), ' uA'])
disp(['Le courant ICQ de Q3 = ', num2str(Ic3_E3*1000), ' mA'])
disp(['Le courant ICQ de Q2 = ', num2str(Ic2_E3), ' A'])
disp(['La tension VCEQ de Q3 = ', num2str(Vce3_E3), ' V'])
disp(['La tension VCEQ de Q2 = ', num2str(Vce2_E3), ' V'])
disp(' ')
%% Calculs PhotoDiode
disp('=====PhotoDiode======')
rSPS = 1070/10000;  % Rayon de la surface sensible de photodiode (en cm)
Sps = pi*(rSPS^2);  % Surface sensible de photodiode (en cm2)
Dspl = 6.1;         % Densité surfacique de puissance lumineuse (En mW/cm2)
PL = Sps*Dspl/1000; % En W
WL =  850;          % longueur d’onde (En nm)
Resp = 550;         % responsivité (En mA/W)
Ipd = Resp*PL;
% impédance interne de 1 GΩ
% Capacité de 15 pF
% Dark current 5nA
disp(['Surface photosensible de la photodiode = ', num2str(Sps), ' cm2'])
disp(['Le photocourant généré par la photodiode = ', num2str(Ipd), ' mA'])
disp(['La puissance lumineuse déposée dans la photodiode = ', num2str(PL*1000), ' mW'])
disp(' ')
%% Calculs Diode
disp('=====Diode======')
% P.162
Eg = 1.424;         % En eV
ni = 2.18e6;        % En cm^−3
er = 13.1*e0;       % Constante diélectrique relative
NA = 10^17;         % Densité de dopage (En cm^−3)
ND = NA;            % Densité de dopage (En cm^−3)
DN = 184;           % Coefficient de diffusion des électrons (en cm^2/s) 
DP = 8.5;           % Coefficient de diffusion des trous (en cm^2/s)
LnCM = 5.09e-5;     % Longueur de diffusion des électrons (En cm)
LpCM = 1.61e-5;     % Longueur de diffusion des trous (En cm)
Ln = LnCM*10000;    % Longueur de diffusion des électrons (En um)
Lp = LpCM*10000;    % Longueur de diffusion des trous (En um)
rAmm = 1.3/2;       % Rayon de la surface active (en mm)
rAcm = rAmm*0.1;    % Rayon de la surface active (en cm)
A = pi*(rAcm^2);    % Dimension de la surface active(en cm^2)
N_Diode = 3.4;      % Facteur d’idéalité (efficacité d’émission)
Rdiode = 0.25;      % Résistance série parasite de (En Ω)

Is_DEL = A*q*(ni^2)*((DP/(LpCM*ND))+(DN/(LnCM*NA)));
V0 = Vt*log((NA*ND)/(ni^2));
Cj0 = A*sqrt(((er*q)/2)*((NA*ND)/(NA+ND))*(1/V0));

disp(['Le courant de saturation de la DEL = ', num2str(Is_DEL), ' A'])
disp(['Le potentiel de contact = ', num2str(V0),' V'])
disp(['La capacité de jonction Cj0 = ', num2str(Cj0*1000000000),' nF'])
disp(' ')