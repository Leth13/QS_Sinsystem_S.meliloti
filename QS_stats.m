
clear all;
close all;
% parameters and rates

numberofcells = 80; 
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
DNAExpR_init = 1;

parameters =  [beta, alpha, deltaI, deltaO, mu, epsilon, k2, k1, k3, k4, k5, k6, k7, epsilon_sinr, rho];
iter = 30;
Initials = [EXPR_init, SINR_init, SINI_init, EXPRstar_init, AHL_iinit, AHL_oinit, DNAsinRstar_init, DNAsinIstar_init, DNAsinR_init, DNAsinI_init, DNAExpR_init];
realpulses = size(iter);
for i = 1:iter
    Initials = [EXPR_init, SINR_init, SINI_init, EXPRstar_init, AHL_iinit, AHL_oinit, DNAsinRstar_init, DNAsinIstar_init, DNAsinR_init, DNAsinI_init, DNAExpR_init];
    [SINIplot, SINRplot, EXPRstarplot, EXPRplot, DNAsinRstarplot, timeplot, meantime] = SSA_QS(Initials, parameters, numberofcells, numberofreactions);
    npeaks = size(numberofcells);
    for t = 1:numberofcells
        [peakLoc] = peakfinder(SINIplot(:,t), ((max(SINIplot(:,t))-min(SINIplot(:,t)))/4), 1);
        npeaks(t) = length(peakLoc);
    end
    realpulses(i) = (mean(npeaks)/(meantime))*60;
end
meanpulses = mean(mean(realpulses));


realpulses_expr = size(iter);
Initials(end) = 0;
for i = 1:iter
    [SINIplot, SINRplot, EXPRstarplot, EXPRplot, DNAsinRstarplot, timeplot, meantime] = SSA_QS(Initials, parameters, numberofcells, numberofreactions);
    npeaks = size(numberofcells);
    for t = 1:numberofcells
        [peakLoc] = peakfinder(SINIplot(:,t), ((max(SINIplot(:,t))-min(SINIplot(:,t)))/4), 1);
        npeaks(t) = length(peakLoc);
    end
    realpulses_expr(i) = (mean(npeaks)/(meantime))*60;
end
meanpulses_expr = mean(mean(realpulses_expr));
d = table([realpulses, realpulses_expr], categorical({'Wildtype', 'Wildtype','Wildtype','Wildtype','Wildtype','Wildtype', 'Wildtype','Wildtype','Wildtype','Wildtype','Wildtype', 'Wildtype','Wildtype','Wildtype','Wildtype','Wildtype', 'Wildtype','Wildtype','Wildtype','Wildtype','Wildtype', 'Wildtype','Wildtype','Wildtype','Wildtype','Wildtype', 'Wildtype','Wildtype','Wildtype','Wildtype', 'Mutant','Mutant','Mutant','Mutant','Mutant','Mutant','Mutant','Mutant','Mutant','Mutant','Mutant','Mutant','Mutant','Mutant','Mutant','Mutant','Mutant','Mutant','Mutant','Mutant','Mutant','Mutant','Mutant','Mutant','Mutant','Mutant','Mutant','Mutant','Mutant','Mutant'}));
figure(1)
set(gca,'Fontsize',18);
boxplot(d.Var1, d.Var2)
title('WT and $EXPR^-$ pulses', 'Interpreter','latex')
xlabel('Genotype')
ylabel('Pulses per hour')