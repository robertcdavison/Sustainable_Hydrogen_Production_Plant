% This script is used in conjunction with the process modeling code in order to create all plots in the report.

clear;
clc;
close all
 
Price = 1300;
 
 
NPVPmax1 = -Inf;
NPVPmax2 = -Inf;
 
T1 = 873:20:1273;
T2 = 473:20:773;
 
for i = 1:length(T1)
    [X11(i),~,~,~,~,~,~,~,~,~,NPV1(i),NPVP1(i)] = PlantSim(1,1,T1(i),773,20,2.5,0.3,Price);
    
    if NPVP1(i) > NPVPmax1
        NPVPmax1 = NPVP1(i);
    end
end
 
 
for i = 1:length(T1)
    [X12(i),~,~,~,~,~,~,~,~,~,NPV2(i),NPVP2(i)] = PlantSim(1,1,T1(i),773,15,2.5,0.3,Price);
    
    if NPVP2(i) > NPVPmax1
        NPVPmax1 = NPVP2(i);
    end
end
 
 
for i = 1:length(T1)
    [X13(i),~,~,~,~,~,~,~,~,~,NPV3(i),NPVP3(i)] = PlantSim(1,1,T1(i),773,10,2.5,0.3,Price);
    
    if NPVP3(i) > NPVPmax1
        NPVPmax1 = NPVP3(i);
    end
end
 
 
for i = 1:length(T2)
    [X14(i),~,~,~,~,~,~,~,~,~,NPV4(i),NPVP4(i)] = PlantSim(1,1,1153,T2(i),20,2.5,0.3,Price);
 
    if NPVP4(i) > NPVPmax2
        NPVPmax2 = NPVP4(i);
    end
end
 
 
for i = 1:length(T2)
    [X15(i),~,~,~,~,~,~,~,~,~,NPV5(i),NPVP5(i)] = PlantSim(1,1,1153,T2(i),15,2.5,0.3,Price);
 
    if NPVP5(i) > NPVPmax2
        NPVPmax2 = NPVP5(i);
    end
end
 
 
for i = 1:length(T2)
    [X16(i),~,~,~,~,~,~,~,~,~,NPV6(i),NPVP6(i)] = PlantSim(1,1,1153,T2(i),10,2.5,0.3,Price);
 
    if NPVP6(i) > NPVPmax2
        NPVPmax2 = NPVP6(i);
    end
end
 
MaxLine1 = NPVPmax1*ones(length(T1));
MaxLine2 = NPVPmax2*ones(length(T2));
 
figure('DefaultAxesFontSize',12)
plot(NPV1/1000000,NPVP1,'LineWidth',2)
hold on
plot(NPV2/1000000,NPVP2,'LineWidth',2)
hold on
plot(NPV3/1000000,NPVP3,'LineWidth',2)
hold on
plot(NPV4/1000000,NPVP4,'LineWidth',2)
hold on
plot(NPV5/1000000,NPVP5,'LineWidth',2)
hold on
plot(NPV6/1000000,NPVP6,'LineWidth',2)
legend('P = 20 bar','P = 15 bar','P = 10 bar')
xlabel('NPV_p_r_o_j, $MM')
ylabel('NPV_%')
xlim([-85 -20])
 
figure('DefaultAxesFontSize',12)
plot(T1-273,NPVP1,'LineWidth',2)
hold on
plot(T1-273,NPVP2,'LineWidth',2)
hold on
plot(T1-273,NPVP3,'LineWidth',2)
hold on
plot(T1-273,MaxLine1,'g--','LineWidth',2)
legend('P = 20 bar','P = 15 bar','P = 10 bar','Max NPV_%')
xlabel('Reactor 1 Temperature (Deg C)')
ylabel('NPV_%')
 
figure('DefaultAxesFontSize',12)
plot(T2-273,NPVP4,'LineWidth',2)
hold on
plot(T2-273,NPVP5,'LineWidth',2)
hold on
plot(T2-273,NPVP6,'LineWidth',2)
hold on
plot(T2-273,MaxLine2,'g--','LineWidth',2)
legend('P = 20 bar','P = 15 bar','P = 10 bar','Max NPV_%')
xlabel('Reactor 2 Temperature (Deg C)')
ylabel('NPV_%')

clear;
clc;
close all
 
Price = 1300;
 
 
NPVPmax1 = -Inf;
NPVPmax2 = -Inf;
 
X_1 = 0.01:0.01:0.99;
X_2 = 0.01:0.01:0.99;
 
for i = 1:length(X_1)
    [X11(i),~,~,~,~,~,~,~,~,~,NPV1(i),NPVP1(i)] = PlantSim(X_1(i),1,1153,573,20,2.5,0.3,Price);
    [X12(i),~,~,~,~,~,~,~,~,~,NPV2(i),NPVP2(i)] = PlantSim(X_1(i),1,1153,573,20,3.0,0.3,Price);
    [X13(i),~,~,~,~,~,~,~,~,~,NPV3(i),NPVP3(i)] = PlantSim(X_1(i),1,1153,573,20,3.5,0.3,Price);
 
    if NPVP1(i) > NPVPmax1
        NPVPmax1 = NPVP1(i);
    end
    if NPVP2(i) > NPVPmax1
        NPVPmax1 = NPVP2(i);
    end
    if NPVP3(i) > NPVPmax1
        NPVPmax1 = NPVP3(i);
    end
end
 
 
for i = 1:length(X_2)
    [~,X21(i),~,~,~,~,~,~,~,~,NPV4(i),NPVP4(i)] = PlantSim(1,X_2(i),1153,573,20,2.5,0.3,Price);
    [~,X22(i),~,~,~,~,~,~,~,~,NPV5(i),NPVP5(i)] = PlantSim(1,X_2(i),1153,573,20,3.0,0.3,Price);
    [~,X33(i),~,~,~,~,~,~,~,~,NPV6(i),NPVP6(i)] = PlantSim(1,X_2(i),1153,573,20,3.5,0.3,Price);
    if NPVP4(i) > NPVPmax2
        NPVPmax2 = NPVP4(i);
    end
    if NPVP5(i) > NPVPmax2
        NPVPmax2 = NPVP5(i);
    end
    if NPVP6(i) > NPVPmax2
        NPVPmax2 = NPVP6(i);
    end
end
 
 
MaxLine1 = NPVPmax1*ones(length(X_1));
MaxLine2 = NPVPmax2*ones(length(X_2));
 
figure('DefaultAxesFontSize',12)
plot(X11,NPVP1,'LineWidth',2)
hold on
plot(X12,NPVP2,'LineWidth',2)
hold on
plot(X13,NPVP3,'LineWidth',2)
hold on
plot(X_1,MaxLine1,'g--','LineWidth',2)
legend('MR = 2.5','MR = 3.0','MR = 3.5','Max NPV_%')
xlabel('Reactor 1 Conversion (X_C_H_4)')
ylabel('NPV_%')
 
figure('DefaultAxesFontSize',12)
plot(X21,NPVP4,'LineWidth',2)
hold on
plot(X22,NPVP5,'LineWidth',2)
hold on
plot(X33,NPVP6,'LineWidth',2)
hold on
plot(X_2,MaxLine2,'g--','LineWidth',2)
legend('MR = 2.5','MR = 3.0','MR = 3.5','Max NPV_%')
xlabel('Reactor 2 Conversion (X_C_O)')
ylabel('NPV_%')

clear;
clc;
close all;
 
Price = 1300;
dW = 0.3;
X_1 = 0.01:0.01:0.99;
X_2 = X_1;
P = 20;
 
for i = 1:length(X_1)
    [X1(i),X2(i),Reactor1Comp(i,:),Reactor2Comp(i,:),Vtot1(i),Vtot2(i),Reactor1Prod(i,:),Reactor2Prod(i,:),FT_1(i),FT_2(i),NPV(i),NPVP(i)] = PlantSim(X_1(i),1,1153,573,P,2.5,dW,Price);
end
 
figure('DefaultAxesFontSize',12)
plot(X1,Reactor1Comp(:,1),X1,Reactor1Comp(:,2),X1,Reactor1Comp(:,3),X1,Reactor1Comp(:,4),X1,Reactor1Comp(:,5),'LineWidth',2)
xlabel('Conversion of Methane (X_C_H_4)')
ylabel('Composition Leaving Reactor 1 (y_i)')
legend('CH_4','H_2O','H_2','CO_2','CO')
 
figure('DefaultAxesFontSize',12)
plot(X1,Vtot1,'LineWidth',2)
xlabel('Conversion of Methane (X_C_H_4)')
ylabel('Reactor 1 Volume (m^3)')
 
figure('DefaultAxesFontSize',12)
plot(X1,Reactor1Prod(:,1),X1,Reactor1Prod(:,2),X1,Reactor1Prod(:,3),X1,Reactor1Prod(:,4),X1,Reactor1Prod(:,5),'LineWidth',2)
xlabel('Conversion of Methane (X_C_H_4)')
ylabel('Production of Reactor 1 (F_i) (mol/s)')
legend('CH_4','H_2O','H_2','CO_2','CO')
 
figure('DefaultAxesFontSize',12)
plot(X1,FT_1,'LineWidth',2)
xlabel('Conversion of Methane (X_C_H_4)')
ylabel('Flowrate into Reactor 1 (mol/s)')
 
figure('DefaultAxesFontSize',12)
plot(X1,Reactor2Comp(:,1),X1,Reactor2Comp(:,2),X1,Reactor2Comp(:,3),X1,Reactor2Comp(:,4),X1,Reactor2Comp(:,5),'LineWidth',2)
xlabel('Conversion of Methane (X_C_H_4)')
ylabel('Composition Into Separation System (y_i)')
legend('CH_4','H_2O','H_2','CO_2','CO')
 
for i = 1:length(X_2)
    [X1(i),X2(i),Reactor1Comp(i,:),Reactor2Comp(i,:),Vtot1(i),Vtot2(i),Reactor1Prod(i,:),Reactor2Prod(i,:),FT_1(i),FT_2(i),NPV(i),NPVP(i)] = PlantSim(1,X_2(i),1153,573,P,2.5,dW,Price);
end
 
figure('DefaultAxesFontSize',12)
plot(X2,Reactor2Comp(:,1),X2,Reactor2Comp(:,2),X2,Reactor2Comp(:,3),X2,Reactor2Comp(:,4),X2,Reactor2Comp(:,5),'LineWidth',2)
xlabel('Conversion of Carbon Monoxide (X_C_O)')
ylabel('Composition Leaving Reactor 2 (y_i)')
legend('CH_4','H_2O','H_2','CO_2','CO')
 
figure('DefaultAxesFontSize',12)
plot(X2,Vtot2,'LineWidth',2)
xlabel('Conversion of Carbon Monoxide (X_C_O)')
ylabel('Reactor 2 Volume (m^3)')
 
figure('DefaultAxesFontSize',12)
plot(X2,Reactor2Prod(:,1),X2,Reactor2Prod(:,2),X2,Reactor2Prod(:,3),X2,Reactor2Prod(:,4),X2,Reactor2Prod(:,5),'LineWidth',2)
xlabel('Conversion of Methane (X_C_H_4)')
ylabel('Production of Reactor 2 (F_i) (mol/s)')
legend('CH_4','H_2O','H_2','CO_2','CO')
 
figure('DefaultAxesFontSize',12)
plot(X2,FT_2,'LineWidth',2)
xlabel('Conversion of Carbon Monoxide (X_C_O)')
ylabel('Flowrate into Reactor 2 (mol/s)')