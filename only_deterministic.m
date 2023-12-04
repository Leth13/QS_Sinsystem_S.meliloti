
clear all;
close all;
% parameters and rates

Parameters_From_File = load("parameters.txt");
beta =    Parameters_From_File(1);% min^-1 - Buceta
alpha =   Parameters_From_File(2); % min^-1 - Production rate of AHL_I - Buceta
deltaI =  Parameters_From_File(3); % min^-1 - infusion rate of AHL_E - Buceta
deltaO =  Parameters_From_File(4); % min^-1 - Diffusion rate of AHL_E - Buceta
mu =      Parameters_From_File(5); %min^-1 - growth rate according to https://doi.org/10.1128/mSphere.00567-18. look at Bettenworth too
epsilon = Parameters_From_File(6); % min^-1  - Protein degradation - Bettenworth
k2 = Parameters_From_File(7); % min^-1 * nM^-1 ExpR - AHL association rate - Buceta
k1 = Parameters_From_File(8); % min^-1 - ExpR* dissociation rate - Buceta
k3 = Parameters_From_File(9); % min^-1 * nM^-1 SINR - DNA_SINI association rate - Buceta
k4 = Parameters_From_File(10); % min^-1 - DNA_sinI* dissociation rate - Buceta
k5 = Parameters_From_File(11); % min^-1 * nM^-1 ExpR* - DNA_SINR association rate - Buceta
k6 = Parameters_From_File(12); % min^-1 - DNA_SINR* dissociation rate - Buceta
k7 = Parameters_From_File(13); % min^-1 * nM^-1 ExpR - SINR - DNA_SINI association rate - Buceta
epsilon_sinr = Parameters_From_File(14); 
rho = Parameters_From_File(15);
volume = 1.7;

%initial quantities of metabolites
EXPR_init = 0;
SINR_init = 0;
SINI_init = 0;
EXPRstar_init = 0;
AHL_iinit = 0;
AHL_oinit = 0;
DNAsinR_init = 1;
DNAsinRstar_init = 0;
DNAsinI_init = 1;
DNAsinIstar_init = 0;


parameters =  [beta, alpha, deltaI, deltaO, mu, epsilon, k2, k1, k3, k4, k5, k6, k7, epsilon_sinr, rho];

Initials = [EXPR_init, SINR_init, SINI_init, EXPRstar_init, AHL_iinit, AHL_oinit, DNAsinRstar_init, DNAsinIstar_init, DNAsinR_init, DNAsinI_init];
temp = [0 2500];
[tdet, det] = ode45(@(t,y) myode_full(y , parameters), temp, Initials);

xxis = 2500;



subplot(3, 2, 1)
plot(tdet, det(:,1));
title("EXPR")
xlabel('time [Min]','interpreter','latex');
ylabel('concentrarion [nmol]','interpreter','latex')
yxis = max(max(det(:,1)))*1.007;
axis([0 xxis 0 yxis])

subplot(3, 2, 2)
plot(tdet, det(:,2));
title("SINR")
xlabel('time [Min]','interpreter','latex');
ylabel('concentrarion [nmol]','interpreter','latex');
yxis = max(max(det(:,2)))*1.5;
axis([0 500 0 yxis])


subplot(3, 2, 3)
plot(tdet, det(:,3));
title("SINI")
xlabel('time [Min]','interpreter','latex');
ylabel('concentrarion [nmol]','interpreter','latex');
yxis = max(max(det(:,3)))*1.007;
axis([0 xxis 0 yxis])


subplot(3, 2, 4)
plot(tdet, det(:,4));
title("EXPR*")
xlabel('time [Min]','interpreter','latex');
ylabel('concentrarion [nmol]','interpreter','latex');
yxis = max(max(det(:,4)))*1.007;
axis([0 xxis 0 yxis])

subplot(3, 2, 5)
plot(tdet, det(:,5))
title("AHLi")
xlabel('time [Min]','interpreter','latex');
ylabel('concentrarion [nmol]','interpreter','latex');
yxis = max(max(det(:,5)))*1.007;
axis([0 xxis 0 yxis])

subplot(3, 2, 6)
plot(tdet, det(:,6))
title("AHLo")
xlabel('time [Min]','interpreter','latex');
ylabel('concentrarion [nmol]','interpreter','latex');
yxis = max(max(det(:,6)))*1.007;
axis([0 xxis 0 yxis])

%expr- plot
[tdet, det] = ode45(@(t,y) myode_expr(y , parameters), temp, Initials);

figure(2)
subplot(2, 2, 1)
plot(tdet, det(:,2));
title("SINR")
xlabel('time [Min]','interpreter','latex');
ylabel('concentrarion [nmol]','interpreter','latex');
yxis = max(max(det(:,2)))*1.5;
axis([0 500 0 yxis])

subplot(2, 2, 2)
plot(tdet, det(:,3));
title("SINI")
xlabel('time [Min]','interpreter','latex');
ylabel('concentrarion [nmol]','interpreter','latex');
yxis = max(max(det(:,3)))*1.10;
axis([0 xxis 0 yxis])

subplot(2, 2, 3)
plot(tdet, det(:,5))
title("AHLi")
xlabel('time [Min]','interpreter','latex');
ylabel('concentrarion [nmol]','interpreter','latex');
yxis = max(max(det(:,5)))*1.1;
axis([0 xxis 0 yxis])
subplot(2, 2, 4)

plot(tdet, det(:,6))
title("AHLo")
xlabel('time [Min]','interpreter','latex');
ylabel('concentrarion [nmol]','interpreter','latex');
yxis = max(max(det(:,6)))*1.1;
axis([0 xxis 0 yxis])
