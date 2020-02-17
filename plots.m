clear;
clc;
 
T1 = 973:100:1273;
P = 10;
for i = 1:length(T1)
    [X1(i),X2(i),Reactor1Comp(i,:),Reactor2Comp(i,:),Vtot1(i),Vtot2(i),Reactor1Prod(i,:),Reactor2Prod(i,:),FT_1(i),FT_2(i),NPV(i),NPVP(i)] = SIMPLE(1,1,T1(i),773,P,3);
end
figure('DefaultAxesFontSize',12)
plot(T1,NPV,'LineWidth',2)
hold on
P = 15;
for i = 1:length(T1)
    [X1(i),X2(i),Reactor1Comp(i,:),Reactor2Comp(i,:),Vtot1(i),Vtot2(i),Reactor1Prod(i,:),Reactor2Prod(i,:),FT_1(i),FT_2(i),NPV(i),NPVP(i)] = SIMPLE(1,1,T1(i),773,P,3);
end
plot(T1,NPV,'LineWidth',2)
hold on
P = 20;
for i = 1:length(T1)
    [X1(i),X2(i),Reactor1Comp(i,:),Reactor2Comp(i,:),Vtot1(i),Vtot2(i),Reactor1Prod(i,:),Reactor2Prod(i,:),FT_1(i),FT_2(i),NPV(i),NPVP(i)] = SIMPLE(1,1,T1(i),773,P,3);
end
plot(T1,NPV,'LineWidth',2)
xlabel('Temperature in PBR 1 (K)')
ylabel('Net Present Value After Year 12 ($)')
legend('P = 10 bar', 'P = 15 bar', 'P = 20 bar')
 
figure('DefaultAxesFontSize',12)
plot(X1,NPVP,'LineWidth',2)
xlabel('Conversion of Methane (X_C_H_4)')
ylabel('NPV Percent After Year 12 (%)')
 
figure('DefaultAxesFontSize',12)
plot(X1,Reactor1Comp(:,1),X1,Reactor1Comp(:,2),X1,Reactor1Comp(:,3),X1,Reactor1Comp(:,4),X1,Reactor1Comp(:,5),'LineWidth',2)
xlabel('Conversion of Methane (X_C_H_4)')
ylabel('Composition Leaving Reactor 1 (y_i)')
legend('CH_4','H_2O','H_2','CO_2','CO')
 
figure('DefaultAxesFontSize',12)
plot(X2,Reactor2Comp(:,1),X2,Reactor2Comp(:,2),X2,Reactor2Comp(:,3),X2,Reactor2Comp(:,4),X2,Reactor2Comp(:,5),'LineWidth',2)
xlabel('Conversion of Carbon Monoxide (X_C_O)')
ylabel('Composition Leaving Reactor 2 (y_i)')
legend('CH_4','H_2O','H_2','CO_2','CO')
 
figure('DefaultAxesFontSize',12)
plot(X1,Vtot1,'LineWidth',2)
xlabel('Conversion of Methane (X_C_H_4)')
ylabel('Reactor 1 Volume (m^3)')
 
figure('DefaultAxesFontSize',12)
plot(X2,Vtot2,'LineWidth',2)
xlabel('Conversion of Carbon Monoxide (X_C_O)')
ylabel('Reactor 2 Volume (m^3)')
 
figure('DefaultAxesFontSize',12)
plot(X1,Reactor1Prod(:,1),X1,Reactor1Prod(:,2),X1,Reactor1Prod(:,3),X1,Reactor1Prod(:,4),X1,Reactor1Prod(:,5),'LineWidth',2)
xlabel('Conversion of Methane (X_C_H_4)')
ylabel('Production of Reactor 1 (F_i) (mol/s)')
legend('CH_4','H_2O','H_2','CO_2','CO')
 
figure('DefaultAxesFontSize',12)
plot(X2,Reactor2Prod(:,1),X2,Reactor2Prod(:,2),X2,Reactor2Prod(:,3),X2,Reactor2Prod(:,4),X2,Reactor2Prod(:,5),'LineWidth',2)
xlabel('Conversion of Methane (X_C_H_4)')
ylabel('Production of Reactor 2 (F_i) (mol/s)')
legend('CH_4','H_2O','H_2','CO_2','CO')
 
figure('DefaultAxesFontSize',12)
plot(X1,FT_1,'LineWidth',2)
xlabel('Conversion of Methane (X_C_H_4)')
ylabel('Flowrate into Reactor 1 (mol/s)')
 
figure('DefaultAxesFontSize',12)
plot(X2,FT_2,'LineWidth',2)
xlabel('Conversion of Carbon Monoxide (X_C_O)')
ylabel('Flowrate into Reactor 2 (mol/s)')
 
figure('DefaultAxesFontSize',12)
plot(X1,Reactor2Comp(:,1),X1,Reactor2Comp(:,2),X1,Reactor2Comp(:,3),X1,Reactor2Comp(:,4),X1,Reactor2Comp(:,5),'LineWidth',2)
xlabel('Conversion of Methane (X_C_H_4)')
ylabel('Composition Into Separation System (y_i)')
legend('CH_4','H_2O','H_2','CO_2','CO')