%%
clear variables; close all; clc;
load('../../../CORE_SIM_1.2_edited/output/RESULTS.mat')
load('../../../CORE_SIM_1.2_edited/input/XS_data.mat')
load('../../../CORE_SIM_1.2_edited/input/GEOM_data.mat')
load("C:/Users/kriped/Chalmers/XEROM/Matlab code/ROM_Analytical_SC/input/ROM_input.mat","height","radius","power","sigmaX","gammaI", "gammaX", "lambdaI", "lambdaX","kappa");
load("../data/Nodal_Values.mat","XS")
%% create sparse CR matrix

CR_SA1 = zeros(size(D1));
CR_SA2 = zeros(size(D1));
%del_SA1 = 8.9E-4; % delta Sigma_a for the fast group
%del_SA2 = 3.4E-4; % delta Sigma_a for the thermal group

del_SA1 = 8.9; % delta Sigma_a for the fast group
del_SA2 = 3.4; % delta Sigma_a for the thermal group


CR_SA1(12:13,12:13,1:11) = del_SA1; CR_SA2(12:13,12:13,1:11) = del_SA2; % 1st control rod thermal and fast
CR_SA1(12:13,20:21,1:11) = del_SA1; CR_SA2(12:13,20:21,1:11) = del_SA2; % 2nd control rod thermal and fast
CR_SA1(20:21,12:13,1:11) = del_SA1; CR_SA2(20:21,12:13,1:11) = del_SA2; % 3rd control rod thermal and fast
CR_SA1(20:21,20:21,1:11) = del_SA1; CR_SA2(20:21,20:21,1:11) = del_SA2; % 4th control rod thermal and fast
CR_SA1(16:17,28:29,1:11) = del_SA1; CR_SA2(16:17,28:29,1:11) = del_SA2; % 5th control rod thermal and fast
CR_SA1(28:29,16:17,1:11) = del_SA1; CR_SA2(28:29,16:17,1:11) = del_SA2; % 6th control rod thermal and fast
CR_SA1(4:5,16:17,1:11) = del_SA1; CR_SA2(4:5,16:17,1:11) = del_SA2; % 7th control rod thermal and fast
CR_SA1(16:17,4:5,1:11) = del_SA1; CR_SA2(16:17,4:5,1:11) = del_SA2; % 8th control rod thermal and fast


%% Create vectors and matrices
sizes = size(MOD1);
sizex = sizes(1);
sizey = sizes(2);
sizez = sizes(3);
%M = sizes(4);
M=10;
v1 = 2E9; % fast neutron velocity in cm/s
v2 = 2.2E5; % thermal neutron velocity in cm/s 
KAPPA = 0.2976e-10; % Guessed value of Kappa (source unknown) same as is used in the 1G homogenous model
%FB = 3.08228E-19; % original Feedback coeficient Calculated in mathematica from the one group homogenous model drdp*(Sigma_a + D*B^2)Kappa*Int(Phi_eq,V)*SigmaF
FB= 7.19198E-19; % test feedback

V_inv = [1/v1, 0 ; 0, 1/v2];
K_VALUE = lambda; % K values / eigenvalues of the modes
DV = DX*DY*DZ; % Discreet volume element
KN = XS.KN; % kappa / nu
NU = kappa./KN; % space dependent nu based on guess of a non space dependent kappa (APPROXIMATION)
SIGF1 = NUFIS1./NU; % Fast fission cross section
SIGF2 = NUFIS2./NU; % Thermal fission cross section
ZERO = zeros(size(NUFIS1)); % zero element matching the size of the reactor
ONE = ones(size(NUFIS1)); % unit element matching the size of the reactor
CR_SA = [CR_SA1,ZERO;ZERO,CR_SA2];
F = [NUFIS1, NUFIS2;ZERO,ZERO]; % FIssion matrix
%D = [-D1,ZERO;ZERO,-D2]; % Diffusion coefficient matrix
MOD = [MOD1;MOD2]; % vector of solutions to the forward problem
MOD_adj = [MOD1_adj,MOD2_adj]; % vector of solutions to the adjoint problem 
MOD_EQ = [MOD1(:,:,:,1);MOD2(:,:,:,1)]; % Vector of only the equilibrium neutron flux solution
PS = power*K_VALUE(1)/(kappa*DV*sum(G2_inner_product([SIGF1,SIGF2],MOD_EQ,"vector","vector"),"all"));
MOD_eq_MAT = PS*[MOD1(:,:,:,1),ZERO; ZERO,MOD2(:,:,:,1)]; % Diagonal matrix containing the solutions to the forward problem
MOD_EQ_scaled = PS*MOD_EQ; % Vector of only the equilibrium neutron flux solution scaled by the power
MOD_UPPER = PS*[MOD2(:,:,:,1), ZERO ; ZERO, ZERO]; % Costom matrix used in the equations 
MOD_LOWER = PS*[ZERO, ZERO; MOD2(:,:,:,1), ZERO]; % Costom matrix used in the equations

PHI_bottom = zeros(1,M);
PHI_top = zeros(1,M);
for mode = 1:M
    PHI_bottom(mode) = DV*sum(G2_inner_product([SIGF1(:,:,1:17),SIGF2(:,:,1:17)],MOD(:,:,1:17,mode),"vector","vector"),1:3);
    PHI_top(mode) = DV*sum(G2_inner_product([SIGF1(:,:,17+1:end),SIGF2(:,:,17+1:end)], MOD(:,:,17+1:end,mode),"vector","vector"),1:3);
end

%SIG = [ABS1+REM, ZERO ; -REM, ABS2]; % Matrix containing the absorbtion and removal cross sections
GAMMAI = [gammaI*SIGF1,gammaI*SIGF2 ; ZERO, ZERO ]; % matrix containing the production of Iodine from fission
GAMMAX = [gammaX*SIGF1,gammaX*SIGF2 ; ZERO,ZERO ]; % matrix containing the production of xenon from fission
%VECX = [ONE;ZERO]; % transformation vector
VECXT = [ONE,ZERO]; % transformation vector
%MATX = [ZERO,ONE;ZERO,ZERO]; %transformation matrix


I0 = 1/(K_VALUE(1)*lambdaI)*G2_inner_product(VECXT, G2_inner_product(GAMMAI,MOD_EQ_scaled,"matrix","vector"),"vector","vector"); % I_0(r) = \hat(X)^T \cdot 1/(keff*lambda_I)Gamma_I\times Phi_0
X0 = (lambdaI*I0 + G2_inner_product(VECXT,1/K_VALUE(1)*G2_inner_product(GAMMAX,MOD_EQ_scaled,"matrix","vector"),"vector","vector"))./(lambdaX+sigmaX*MOD2(:,:,:,1)); 
Xe_UPPER = [ZERO,X0;ZERO,ZERO]; % Custom matrix used in the equations


clear NUFIS1 NUFIS2 n_iter n_restart ABS1 ABS2 FLX1 FLX2 EIG_MET D1 D2 MOD1_adj MOD2_adj REM RES_FLX RES_MOD RES_MOD_adj ...
    XS KN DX DY DZ conv_ERAM conv_POW lambda lambda_adj gammaI gammaX v1 v2
%% initialise parameters
PHID_F_PHI = zeros(1,M);
PHID_V_PHI = zeros(1,M);
PHID_GAMMAI_PHI = zeros(M);
PHID_GAMMAX_PHI = zeros(M);
PHID_PHI_eq_mat_PHI = zeros(M);
PHID_PHI_eq_PHI = zeros(M);
PHID_PHIUPPER_PHI=zeros(M);
PHID_PHILOWER_PHI=zeros(M);
PHID_X0_PHI=zeros(M);
PHID_CR_PHI = zeros(M);
temp_GAMMAX_PHI = zeros(sizex*2,sizey,sizez,M);
temp_GAMMAI_PHI = zeros(sizex*2,sizey,sizez,M);
temp_F_PHI = zeros(sizex*2,sizey,sizez,M);
temp_V_PHI = zeros(sizex*2,sizey,sizez,M);
temp_PHI_eq_mat_PHI=zeros(sizex*2,sizey,sizez,M);
temp_PHIUPPER_PHI=zeros(sizex*2,sizey,sizez,M);
temp_PHILOWER_PHI=zeros(sizex*2,sizey,sizez,M);
temp_X0_PHI = zeros(sizex*2,sizey,sizez,M);
temp_CR_PHI = zeros(sizex*2,sizey,sizez,M);
%% Calculate kets
for n = 1:M
    temp_F_PHI(:,:,:,n)=G2_inner_product(F(:,:,:),MOD(:,:,:,n),"matrix","vector"); % | F* Phi_n>
    temp_PHI_eq_mat_PHI(:,:,:,n) = G2_inner_product(MOD_eq_MAT,MOD(:,:,:,n),"matrix","vector"); %|Phi_0_mat * Phi_n>
    temp_GAMMAX_PHI(:,:,:,n) = 1/K_VALUE(1)*G2_inner_product(GAMMAX,MOD(:,:,:,n),"matrix","vector"); %  |1/k_0 *Gamma_X * Phi_n>
    temp_GAMMAI_PHI(:,:,:,n) = 1/K_VALUE(1)*G2_inner_product(GAMMAI,MOD(:,:,:,n),"matrix","vector"); % | 1/k_0 Gamma_I * Phi_n>
    temp_PHIUPPER_PHI(:,:,:,n) = G2_inner_product(MOD_UPPER,temp_F_PHI(:,:,:,n),"matrix","vector"); % | \bar{X} \Phi_0 \hat{X}^T * F* Phi_n >
    temp_PHILOWER_PHI(:,:,:,n) = G2_inner_product(MOD_LOWER,temp_F_PHI(:,:,:,n),"matrix","vector"); % | \tilde{X} \Phi_0 \hat{X}^T * F* Phi_n >
    temp_X0_PHI(:,:,:,n) = G2_inner_product(Xe_UPPER,MOD(:,:,:,n),"matrix","vector"); % |\bar{X} * X0 * Phi_n >
    temp_CR_PHI(:,:,:,n) = G2_inner_product(CR_SA,MOD(:,:,:,n),"matrix","vector"); % CP
end
%% Calculate bra-kets 
for m = 1:M
    PHID_F_PHI(m) = DV*sum(G2_inner_product(MOD_adj(:,:,:,m),temp_F_PHI(:,:,:,m),"vector","vector"),'all'); %<Phi^dagger_m| F Phi_m>
    temp_V_PHI(:,:,:,m) = G2_inner_product(V_inv,MOD(:,:,:,m),"scalar_matrix","vector"); % |v^-1 Phi_m>
    PHID_V_PHI(m) = DV*sum(G2_inner_product(MOD_adj(:,:,:,m),temp_V_PHI(:,:,:,m),"vector","vector"),"all"); %<Phi^dagger_m| v^-1 Phi_m>
    for n = 1:M
        PHID_GAMMAX_PHI(m,n) = DV*sum(G2_inner_product(MOD_adj(:,:,:,m),temp_GAMMAX_PHI(:,:,:,n),"vector","vector"),"all"); % <Phi^dagger_m | Gamma_X * Phi_n >
        PHID_GAMMAI_PHI(m,n) =  DV*sum(G2_inner_product(MOD_adj(:,:,:,m),temp_GAMMAI_PHI(:,:,:,n),"vector","vector"),"all"); %<Phi^dagger_m | Gamma_I * Phi_n > 
        PHID_PHI_eq_mat_PHI(m,n) =  DV*sum(G2_inner_product(MOD_adj(:,:,:,m),temp_PHI_eq_mat_PHI(:,:,:,n),"vector","vector"),"all"); %<Phi^dagger_m | Phi_eq_mat Phi_n >
        PHID_PHIUPPER_PHI(m,n) = DV*sum(G2_inner_product(MOD_adj(:,:,:,m),temp_PHIUPPER_PHI(:,:,:,n),"vector","vector"),"all"); %<Phi^dagger_m | \bar{X} \Phi_0 \hat{X}^T * F Phi_n >
        PHID_PHILOWER_PHI(m,n) = DV*sum(G2_inner_product(MOD_adj(:,:,:,m),temp_PHILOWER_PHI(:,:,:,n),"vector","vector"),"all"); %<Phi^dagger_m | \tilde{X} \Phi_0 \hat{X}^T * F Phi_n >
        PHID_X0_PHI(m,n) = DV*sum(G2_inner_product(MOD_adj(:,:,:,m),temp_X0_PHI(:,:,:,n),"vector","vector"),"all"); %<Phi^dagger_m | \bar{X} X_0 Phi_n >
        PHID_CR_PHI(m,n) = DV*sum(G2_inner_product(MOD_adj(:,:,:,m),temp_CR_PHI(:,:,:,n),"vector","vector"),"all"); % <Phi^dagger_m | CR Phi_n >
    end
end
 LAMBDA = PHID_V_PHI./ PHID_F_PHI; % <Phi^dagger_m |v^-1 Phi_n>/ <Phi^dagger_m |F Phi_m> 

clear temp* VECX GAMMAI GAMMAX F D ONE ZERO I0 X0 SIGF1 SIGF2 SIG Xe_UPPER MATX NU

% %% simulate time before perturbation
% 
% tc = 2; % time (in hours) of control rod insertion
% 
% %f = FunctionGen_w_CR(M);
% %fhandle = eval(['@(t,s)[', f, ']']);
% ti=0; tf = tc*3600;
% tspan = [ti,tf];
% %ExMode = 2;
% %exmodeidx = (ExMode-1)*3+1; % index of the exited flux mode
% IC = zeros(1,3*M);
% %IC(exmodeidx) = 0.01;
% opts=odeset("MaxStep",10); % Options for solver
% 
% [time_b, state_values_b] = ode15s(@ode,tspan,IC,opts);
% 

%% evaluate function after introduction of of control rods

% 
% fhandle = eval(['@(t,s)[', f, ']']);


% ti=0; tf = 30*3600;
% tspan = [ti,tf];
% ExMode = 2;
% exmodeidx = (ExMode-1)*3+1; % index of the exited flux mode
% IC = zeros(1,3*M);
% IC(exmodeidx) = 0.01;
% opts=odeset("MaxStep",10); % Options for solver
% 
% [time, state_values] = ode15s(ode,tspan,IC,opts,M);
clear MOD_* CR* DV size*
save ../data/tempFile
[time,state_values_2G]=ivp(M);
save("../data/amplitudes.mat","time","state_values_2G","FB")
%% calculate ASI

% time_total = [time_b;time_a];
% state_values_total = [state_values_a;state_values_b];
% 
% ASI = ((PHI_bottom(1)-PHI_top(1))+sum(state_values_total(:,1:3:end).*(PHI_bottom-PHI_top)))./((PHI_bottom(1)+PHI_top(1))+sum(state_values_total(:,1:3:end).*(PHI_bottom+PHI_top)));  
% 
% figure()

%plot(time, state_values(:,((1:10)-1)*3+1))
% %% Plot flux iodine and xenon
% figure(1)
% timehours = time/3600;
% plotmode = 1;
% plot(timehours, state_values(:,(plotmode-1)*3+1:(plotmode-1)*3+3))
% yl = ylabel("Amplitude [AU]",'Fontsize', 14,"Rotation",90);
% %yl.Position(1) = 78; 
% legend(["flux","Iodine","xenon"],"Location", "northeastoutside")
% %legend(["Mode 2" "Mode 15" "Mode 1" "Mode 7" "Mode 12"],"Location","North")
% xlabel("Time [h]",'Fontsize', 14)
% 
% %ylim([0 1E-6]);
% grid on
% hold off
% 
% 
% %% Plot times signal
% 
% sol_phi = state_values(((1:10)-1)*3+1);
% MOD1_point(1:10) = reshape(MOD1(22,22,22,:),[1,10]);
% %MOD2_point(1:10) = MOD2(22,22,22,:);
% %Time_signal1 = sol_phi .* MOD1_point;
% %Time_signal2 = sol_phi .* MOD2_point;
% 
% %figure(2)
% 
% %plot(time, Time_signal1(2))
% 
% function DY = ode(t,s,M)
%     f = FunctionGen(M);
%     %CR = zeros(M).*(t<=tc) + PHID_CR_PHI.*(t>tc);
%     DY = eval(['[',f, ']']);
% end
