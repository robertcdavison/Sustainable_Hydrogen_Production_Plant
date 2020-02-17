% The following script was used to iterate through all possible design variables in order to find the optimum NPV%

clear;
clc;
close all
dW = 0.3;
Price = 1300;
T1 = 953:20:1273;
T2 = 473:20:773;
X_1 = 0.1:0.01:0.99;
X_2 = 0.1:0.01:0.99;
P = 5:0.5:20;
MR = 2.5:0.1:5;
NPVmax = -Inf;
for i = 1:length(T1)
    [X1,X2,Reactor1Comp,Reactor2Comp,Vtot1,Vtot2,Reactor1Prod,Reactor2Prod,FT_1,FT_2,NPV,NPVP] = PlantSim(1,1,T1(i),773,20,2.5,dW,Price);
    if NPV >= NPVmax
        OptT1 = T1(i);
        NPVmax = NPV;
    end
end
for i = 1:length(T2)
    [X1,X2,Reactor1Comp,Reactor2Comp,Vtot1,Vtot2,Reactor1Prod,Reactor2Prod,FT_1,FT_2,NPV,NPVP] = PlantSim(1,1,OptT1,T2(i),20,2.5,dW,Price);
    if NPV >= NPVmax
        OptT2 = T2(i);
        NPVmax = NPV;
    end
end
for i = 1:length(P)
    [X1,X2,Reactor1Comp,Reactor2Comp,Vtot1,Vtot2,Reactor1Prod,Reactor2Prod,FT_1,FT_2,NPV,NPVP] = PlantSim(1,1,OptT1,OptT2,P(i),2.5,dW,Price);
    if NPV >= NPVmax
        OptP = P(i);
        NPVmax = NPV;
    end
end
for i = 1:length(MR)
    [X1,X2,Reactor1Comp,Reactor2Comp,Vtot1,Vtot2,Reactor1Prod,Reactor2Prod,FT_1,FT_2,NPV,NPVP] = PlantSim(1,1,OptT1,OptT2,OptP,MR(i),dW,Price);
    if NPV >= NPVmax
        OptMR = MR(i);
        NPVmax = NPV;
    end
end
for i = 1:length(X_1)
    [X1,X2,Reactor1Comp,Reactor2Comp,Vtot1,Vtot2,Reactor1Prod,Reactor2Prod,FT_1,FT_2,NPV,NPVP] = PlantSim(X_1(i),1,OptT1,OptT2,OptP,OptMR,dW,Price);
    if NPV >= NPVmax
        OptX1 = X1;
        NPVmax = NPV;
    end
end
for i = 1:length(X_2)
    [X1,X2,Reactor1Comp,Reactor2Comp,Vtot1,Vtot2,Reactor1Prod,Reactor2Prod,FT_1,FT_2,NPV,NPVP] = PlantSim(OptX1,X_2(i),OptT1,OptT2,OptP,OptMR,dW,Price);
    if NPV >= NPVmax
        OptX2 = X2;
        NPVmax = NPV;
    end
end