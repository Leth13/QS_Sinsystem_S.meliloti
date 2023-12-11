function [SINIplot, SINRplot, EXPRstarplot, EXPRplot, timeplot, meantime] = SSA_QS(y, parameters, numberofcells, numberofreactions)

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
volume = 1.2;

EXPR_init = y(1);
SINR_init = y(2);
SINI_init = y(3);
EXPRstar_init = y(4);
AHL_iinit = y(5);
AHL_oinit = y(6);
DNAsinRstar_init = y(7);
DNAsinIstar_init = y(8);
DNAsinR_init = y(9);
DNAsinI_init = y(10);
DNAexpR_init = y(11);

SINRplot = zeros(numberofreactions, numberofcells);
SINIplot = zeros(numberofreactions, numberofcells);
EXPRstarplot = zeros(numberofreactions, numberofcells);
DNAsinRstarplot = zeros(numberofreactions, numberofcells);
EXPRplot = zeros(numberofreactions, numberofcells);
timeplot = zeros(numberofreactions, numberofcells);
totaltime = 0;
for i=1:numberofcells
    EXPR = EXPR_init;
    SINR = SINR_init;
    SINI = SINI_init;
    AHL_i = AHL_iinit;
    AHL_o = AHL_oinit;
    EXPRstar = EXPRstar_init;
    DNAsinR = DNAsinR_init;
    DNAsinRstar = DNAsinRstar_init;
    DNAsinI = DNAsinI_init;
    DNAsinIstar = DNAsinIstar_init;
    DNAexpR = DNAexpR_init;
    step = 0;
    k = 1;
    while (k<=numberofreactions) 
                rr = rand(2,1);
                a1 = beta*DNAexpR; %can't be 0
                a2 = mu*EXPR; 
                a3 = epsilon*EXPR; 
                a4 = beta*DNAsinIstar; 
                a5 = mu*SINI; 
                a6 = epsilon*SINI; 
                a7 = alpha*SINI; 
                a8 = mu*AHL_i;  
                a9 = epsilon*AHL_i; 
                a10 = deltaO*AHL_i; 
                a11 = deltaI*AHL_o; 
                a12 = (k2/volume)*AHL_i*EXPR; 
                a13 = k1*EXPRstar; 
                a14 = mu*EXPRstar; 
                a15 = (k3/volume)*DNAsinI*SINR; 
                a16 = k4*DNAsinIstar;
                a17 = ((k7/volume)/volume)*DNAsinI*SINR*EXPRstar; 
                a18 = k5/volume*EXPRstar*DNAsinR; 
                a19 = k6*DNAsinRstar; 
                a20 = beta*DNAsinR; % can't be 0 at the start
                a21 = epsilon_sinr*SINR;
                a22 = mu*SINR;
                a23 = rho*AHL_o;
                a24 = EXPRstar*epsilon;
                a0 = a1 + a2 + a3 +a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12 + a13 + a14 + a15 + a16 + a17 + a18 + a19 + a20 + a21 + a22 + a23 + a24; 
                step = step+(1/a0)*log(1/rr(1));
                m = rr(2)*a0;
                if m <= a1 
                %beta*DNAexpR
                    EXPR = EXPR + 1;
                elseif m <= (a1 + a2)
                %mu*EXPR
                    EXPR = EXPR - 1;
                elseif m <= (a1 + a2 + a3)
                %epsilon*EXPR
                    EXPR = EXPR -1;
                elseif m <= (a1 + a2 + a3 + a4)
                %beta*DNAsinIstar
                    SINI = SINI + 1;
                elseif m <= (a1 + a2 + a3 + a4 + a5)
                %mu*SINI
                    SINI = SINI - 1;
                elseif m <= (a1 + a2 + a3 + a4 + a5 + a6)
                 %epsilon*SINI
                    SINI = SINI - 1;
                elseif m <= (a1 + a2 + a3 + a4 + a5 + a6 + a7)
                %alpha*SINI
                    AHL_i = AHL_i + 1;
                elseif m <= (a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8)
                 %mu*AHL_i
                    AHL_i = AHL_i - 1;
                elseif m <= (a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9)
                %epsilon*AHL_i
                    AHL_i = AHL_i - 1;
                elseif m <= (a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10)
                %deltaO*AHL_i
                    AHL_i = AHL_i -1;
                    AHL_o = AHL_o +1;
                elseif m <= (a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11)
                %deltaI*AHL_o
                    AHL_i = AHL_i +1;
                    AHL_o = AHL_o -1;
                elseif m < (a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12)
                %k2/volume*AHL_i*EXPR
                    AHL_i = AHL_i - 1;
                    EXPR = EXPR - 1;
                    EXPRstar = EXPRstar + 1;
                elseif m <= (a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12 + a13)
                %k1*EXPRstar
                    AHL_i = AHL_i + 1;
                    EXPR = EXPR + 1;
                    EXPRstar = EXPRstar - 1;
                elseif m <= (a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12 + a13 + a14)
                %mu*EXPRstar
                    EXPRstar = EXPRstar - 1;
                elseif m <= (a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12 + a13 + a14 + a15)
                %k3/volume*DNAsinI*SINR
                    DNAsinIstar = DNAsinIstar + 1;
                    DNAsinI = DNAsinI - 1;
                    SINR = SINR - 1;
                elseif m <= (a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12 + a13 + a14 + a15 + a16)
                %k4*DNAsinIstar
                    DNAsinIstar = DNAsinIstar - 1;
                    DNAsinI = DNAsinI + 1;
                    SINR = SINR + 1;
                elseif m <= (a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12 + a13 + a14 + a15 + a16 + a17)
                %k7/volume/volume*DNAsinI*SINR*EXPRstar
                    DNAsinIstar = DNAsinIstar + 1;
                    DNAsinI = DNAsinI - 1;
                    SINR = SINR - 1;
                elseif m <= (a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12 + a13 + a14 + a15 + a16 + a17 + a18)
                %k5/volume*EXPRstar*DNAsinR
                    DNAsinRstar = DNAsinRstar + 1;
                    DNAsinR = DNAsinR - 1;
                    EXPRstar = EXPRstar - 1;
                elseif m <= (a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12 + a13 + a14 + a15 + a16 + a17 + a18 + a19)
                %k6*DNAsinRstar
                    DNAsinRstar = DNAsinRstar - 1;
                    DNAsinR = DNAsinR + 1;
                    EXPRstar = EXPRstar + 1;
                elseif m <= (a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12 + a13 + a14 + a15 + a16 + a17 + a18 + a19 + a21)
                %beta*DNAsinR%epsilon*SINR
                    SINR = SINR - 1;
                elseif m <= (a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12 + a13 + a14 + a15 + a16 + a17 + a18 + a19 + a21 + a22)
                %mu*SINR
                    SINR = SINR - 1;
                elseif m <= (a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12 + a13 + a14 + a15 + a16 + a17 + a18 + a19 + a21 + a22 + a20)
                %beta*DNAsinR
                    SINR = SINR + 1;
                elseif m <= (a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12 + a13 + a14 + a15 + a16 + a17 + a18 + a19 + a21 + a22 + a20 + a23)
                %rho*AHL_o
                    AHL_o = AHL_o - 1;
                else
                %EXPRstar*epsilon
                    EXPRstar = EXPRstar - 1;
                end
                SINIplot(k, i) = SINI;
                EXPRplot(k, i) = EXPR;
                SINRplot(k, i) = SINR;
                EXPRstarplot(k, i) = EXPRstar;
                timeplot(k, i) = step;
                k = k + 1;
    end
    totaltime = max(timeplot(:, i)) + totaltime;
end

%calculate the mean time of all the iterations
meantime = totaltime/numberofcells;

end
