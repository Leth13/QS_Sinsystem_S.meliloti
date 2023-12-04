function dydt = myode_full(y, parameters)

beta = parameters(1);
alpha = parameters(2);
deltaI = parameters(3);
deltaO = parameters(4);
mu = parameters(5);
epsilon = parameters(6);
k2 = parameters(7);
k1 = parameters(8);
k3 = parameters(9);
k4 = parameters(10);
k5 = parameters(11);
k6 = parameters(12);
k7 = parameters(13);
epsilon_sinr = parameters(14);
rho = parameters(15);
EXPR = y(1);
SINR = y(2);
SINI = y(3);
EXPRstar = y(4);
AHL_i = y(5);
AHL_o = y(6);
DNAsinRstar = y(7);
DNAsinIstar = y(8);
DNAsinI = y(9);
DNAsinR = y(10);
DNAexpR = 1;


dEXPRdt = -k2*EXPR*AHL_i + k1*EXPRstar + beta*DNAexpR -mu*EXPR - epsilon*EXPR;

dEXPRstardt = k2*EXPR*AHL_i - k1*EXPRstar -k5*EXPRstar*DNAsinR + k6*DNAsinRstar - mu*EXPRstar -epsilon*EXPRstar;

dSINRdt = beta*DNAsinR -k3*DNAsinI*SINR + k4*DNAsinIstar - k7*DNAsinI*EXPRstar*SINR - mu*SINR -epsilon_sinr*SINR;

dSINIdt = beta*DNAsinIstar - mu*SINI - epsilon*SINI;

dAHL_idt = alpha*SINI - mu*AHL_i - epsilon*AHL_i - deltaO*AHL_i + deltaI*AHL_o -k2*EXPR*AHL_i + k1*EXPRstar;

dAHL_odt = deltaO*AHL_i - deltaI*AHL_o -  rho*AHL_o;

DnasinRstardt = k5*DNAsinR*EXPRstar - k6*DNAsinRstar;

DNAsinRdt = -k5*DNAsinR*EXPRstar + k6*DNAsinRstar;

DnasinIstardt = k3*DNAsinI*SINR - k4*DNAsinIstar + k7*DNAsinI*EXPRstar*SINR;

DNAsinIdt = -k3*DNAsinI*SINR + k4*DNAsinIstar - k7*DNAsinI*EXPRstar*SINR;


dydt = [dEXPRdt; dSINRdt; dSINIdt; dEXPRstardt;   dAHL_idt; dAHL_odt; DnasinRstardt; DnasinIstardt; DNAsinIdt; DNAsinRdt];
end
