
clear all;
close all;
% parameters and rates

numberofcells = 20; 
numberofreactions = 1000000;


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


%initial quantities of metabolites or number of genes
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
DNAExpR_init = 1;


parameters =  [beta, alpha, deltaI, deltaO, mu, epsilon, k2, k1, k3, k4, k5, k6, k7, epsilon_sinr];

Initials = [EXPR_init, SINR_init, SINI_init, EXPRstar_init, AHL_iinit, AHL_oinit, DNAsinRstar_init, DNAsinIstar_init, DNAsinR_init, DNAsinI_init, DNAExpR_init];
[SINIplot, SINRplot, EXPRplot, timeplot, meantime] = SSA_QS_maxAHL(Initials, parameters, numberofcells, numberofreactions); 

%to start the deterministic model
temp = [0 2500];
% remove DNA_ExpR_init, it's an internal variable of the myode_full
Initials(end) = [];
[tdet, det] = ode45(@(t,y) myode_full(y , parameters), temp, Initials);


%visualizing the models

figure("name", "SINI expression")
set(gca,'Fontsize',18);
for l = 1:numberofcells
    plot(timeplot(:,l),SINIplot(:,l),'HandleVisibility','off');
    hold on
end
hold on;
plot(tdet,det(:,3),'--k','Linewidth',4);
xlabel('time [min]','interpreter','latex');
ylabel('number of $SINI$ molecules','interpreter','latex');
hh=legend('solution of ODEs');
set(hh,'interpreter','latex','location','northwest','Fontsize',18);
yxis = max(max(SINIplot(:,l)))*1.007;
axis([0 xxis 0 yxis])


figure("name", "SINR expression")
set(gca,'Fontsize',18);
for l = 1:numberofcells
    plot(timeplot(:,l),SINRplot(:,l),'HandleVisibility','off');
    hold on
end
hold on;
plot(tdet,det(:,2),'--k','Linewidth',4);
xlabel('time [min]','interpreter','latex');
ylabel('number of $SINR$ molecules','interpreter','latex');
hh=legend('solution of ODEs');
set(hh,'interpreter','latex','location','northwest','Fontsize',18);
yxis = max(max(SINRplot(:,l)))*1.007;
axis([0 xxis 0 yxis])


