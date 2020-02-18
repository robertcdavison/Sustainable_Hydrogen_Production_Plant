% This function will take a desired conversion in reactor 1 (X_1), conversion in reactor 2 (X_2), T1, T2, P, and MR, and output all reactor sizes, an estimate for installation costs, and energy costs.

% It automatically adjusts feed rate so 63 kta of hydrogen is always produced and it will stop increasing reactor volume once a maximum conversion is reached.

% It also calculates NPV based on correlations for equipment sizes and other costs.

function [X1,X2,Reactor1Comp,Reactor2Comp,Vtot1,Vtot2,Reactor1Prod,Reactor2Prod,FT_1,FT_2,NPV,NPVP] = PlantSim(X_1,X_2,T1,T2,P,MR,dW,H2Price)
R1 = 8.314*10^(-3); %R constant for reaction kinetics
format shortEng
format compact
%% Information from Hysys
A1 = 649.3; %Area of each heat exchanger
%% Begin capacity loop
Capacity2 = 1; 
err = 1;
iter = 0;
F1 = 1000;
while err > 0.001
%% Flowrates into first reactor
F1 = F1/Capacity2; %Flow rate of methane feed (mol/s)
F2 = F1*MR; %Inlet feed of H2O (mol/s)
W1 = 0; %Initial catalyst weight in reactor 1 (kg)
W2 = 0; %Initial catalyst weight in reactor 2 (kg)
FCH4 = 0.9963*F1; %molar flow of CH4
FCO2 = 0.0037*F1; %molar flow of CO2
FH2 = 0; %molar flow of H2
FCO = 0; %molar flow of CO
FH2O = F2; %molar flow of H2O
%% Define flow stream into Reactor 1
FCH4_1 = FCH4;
FH2O_1 = FH2O;
FH2_1 = FH2;
FCO_1 = FCO;
FCO2_1 = FCO2;
FT_1 = FCH4 + FH2O + FH2 + FCO2 + FCO;
MdotCH41 = (FCH4*16.04 + FH2*2.016 + FCO2*44.01 + FCO*28)*3600/1000;
MdotH2O1 = (FH2O*18.01)*3600/1000;
MdotR1in = (FCH4*16.04 + FH2*2.016 + FCO2*44.01 + FCO*28 + FH2O*18.01)*3600/1000;
w0CH4 = FCH4*16.04/MdotR1in*3600/1000;
w0H2O = FH2O*18.01/MdotR1in*3600/1000;
w0H2 = FH2*2.016/MdotR1in*3600/1000;
w0CO2 = FCO2*44.01/MdotR1in*3600/1000;
w0CO = FCO*28/MdotR1in*3600/1000;
%% Perform finite difference method to solve Reactor 1
X1 = 0;
while X1 < X_1
    FT = FCH4 + FH2O + FH2 + FCO2 + FCO; % Recalculate total molar flow
    PCH4 = (FCH4/FT)*P; %Recalculate all partial pressures
    PCO = (FCO/FT)*P;
    PH2O = (FH2O/FT)*P;
    PH2 = (FH2/FT)*P;
    PCO2 = (FCO2/FT)*P;
    k1 = k_1(T1,R1); %Calculate rate constants at each point in reactor 1
    k2 = k_2(T1,R1);
    k3 = k_3(T1,R1);
    K1 = K_1(T1);
    K2 = K_2(T1);
    KCO = K_CO(T1,R1); %Calculate equilibrium constants at each point in reactor 1
    KH2 = K_H2(T1,R1);
    KCH4 = K_CH4(T1,R1);
    KH2O = K_H2O(T1,R1);
    PSI = PSI_(KCO,KH2,KCH4,KH2O,PCO,PH2,PCH4,PH2O);
    r1 = r_1(PSI,k1,K1,PH2,PCH4,PH2O,PCO); %Reaction rates at each point in reactor 1
    r2 = r_2(PSI,k2,K2,PH2,PCO2,PH2O,PCO);
    r3 = r_3(PSI,k3,K1,K2,PCH4,PH2,PCO2,PH2O);
    FCH4 = FCH4 - r1*dW - r3*dW; %Change in all component molar flow rates
    FH2O = FH2O - r1*dW - r2*dW - 2*r3*dW;
    FCO = FCO + r1*dW - r2*dW;
    FCO2 = FCO2 + r2*dW + r3*dW;
    FH2 = FH2 + 3*r1*dW + r2*dW + 4*r3*dW;
    X1 = (FCH4_1-FCH4)/(FCH4_1); %Define conversion of methane at each reactor size
    if W1 > 300000
        X_1 = X1;
        break
    end
    W1 = W1 + dW;
end
%% Set up initial conditions flowing into reactor 2
FCH4_2 = FCH4;
FH2O_2 = FH2O;
FH2_2 = FH2;
FCO2_2 = FCO2;
FCO_2 = FCO;
Reactor1Prod = [FCH4_2 FH2O_2 FH2_2 FCO2_2 FCO_2];
FT_2 = FCO_2 + FH2_2 + FCO2_2 + FH2O_2 + FCH4_2;
y1CH4 = FCH4_2/FT_2;
y1H2O = FH2O_2/FT_2;
y1H2 = FH2_2/FT_2;
y1CO2 = FCO2_2/FT_2;
y1CO = FCO_2/FT_2;
Reactor1Comp = [y1CH4 y1H2O y1H2 y1CO2 y1CO];
MdotR1out = (FCH4*16.04 + FH2*2.016 + FCO2*44.01 + FCO*28 + FH2O*18.01)*3600/1000;
w1CH4 = FCH4*16.04/MdotR1out*3600/1000;
w1H2O = FH2O*18.01/MdotR1out*3600/1000;
w1H2 = FH2*2.016/MdotR1out*3600/1000;
w1CO2 = FCO2*44.01/MdotR1out*3600/1000;
w1CO = FCO*28/MdotR1out*3600/1000;
%% Solve Reactor 2 with finite difference
X2 = 0;
while X2 < X_2
    FT = FCH4 + FH2O + FH2 + FCO2 + FCO; % Recalculate total molar flow
    PCO = (FCO/FT)*P; %Recalculate all partial pressures
    PCO2 = (FCO2/FT)*P;
    PH2O = (FH2O/FT)*P;
    PH2 = (FH2/FT)*P;
    k4 = k_4(T2); %Calculate rate constants at each point in reactor 2
    K2 = K_2(T2);
    r4 = r_4(k4,K2,PCO,PH2O,PH2,PCO2); %Reaction rate in reactor 2
    FH2O = FH2O - r4*dW; %Recalculate molar flow rates
    FCO = FCO - r4*dW;
    FCO2 = FCO2 + r4*dW;
    FH2 = FH2 + r4*dW;
    X2 = (FCO_2-FCO)/(FCO_2); %Define conversion of Carbon Monoxide
    if W2 > 300000
        X_2 = X2;
        break
    end
    W2 = W2 + dW; %Add catalyst weight in each iteration
end
MdotR2out = (FCH4*16.04 + FH2*2.016 + FCO2*44.01 + FCO*28 + FH2O*18.01)*3600/1000;
w2CH4 = FCH4*16.04/MdotR2out*3600/1000;
w2H2O = FH2O*18.01/MdotR2out*3600/1000;
w2H2 = FH2*2.016/MdotR2out*3600/1000;
w2CO2 = FCO2*44.01/MdotR2out*3600/1000;
w2CO = FCO*28/MdotR2out*3600/1000;
FCH4_3 = FCH4;
FH2O_3 = FH2O;
FH2_3 = FH2;
FCO2_3 = FCO2;
FCO_3 = FCO;
Reactor2Prod = [FCH4_3 FH2O_3 FH2_3 FCO2_3 FCO_3];
FT_3 = FCO_3 + FH2_3 + FCO2_3 + FH2O_3 + FCH4_3;
y2CH4 = FCH4_3/FT_3;
y2H2O = FH2O_3/FT_3;
y2H2 = FH2_3/FT_3;
y2CO2 = FCO2_3/FT_3;
y2CO = FCO_3/FT_3;
Reactor2Comp = [y2CH4 y2H2O y2H2 y2CO2 y2CO];
FH2O_out = FH2O;
FH2O = 0;
MdotComp = (FCH4*16.04 + FH2*2.016 + FCO2*44.01 + FCO*28 + FH2O*18.01)*3600/1000;
wcompCH4 = FCH4*16.04/MdotComp*3600/1000;
wcompH2O = FH2O*18.01/MdotComp*3600/1000;
wcompH2 = FH2*2.016/MdotComp*3600/1000;
wcompCO2 = FCO2*44.01/MdotComp*3600/1000;
wcompCO = FCO*28/MdotComp*3600/1000;
FT_4 = FH2 + FCH4 + FCO2 + FCO; %Molar flow rate into compressor/PSA (mol/s)
 
%% Solve for size/mass of PSA
Phigh = 30*.987; %High PSA Pressure (atm)
Plow = 2*.987; %Low PSA Pressure (atm)
PCOhigh = Phigh*FCO/FT_4; %Partial pressure of CO into PSA (atm)
PCOlow = Plow*FCO/(FT_4-0.9*FH2); %Partial pressure of CO into PSA (atm)
T = 303; %Temp into PSA
fL = 0.95; %fraction of bed saturation
% Parameters for CO
k1 = 5.05;
k2 = -0.00905;
k3 = 0.001137;
k4 = 1617.0;
k5 = 0.5245;
k6 = 256.5;
n = k5 + k6/T;
B = k3*exp(k4/T);
qm = k1 + k2/T;
qH = (qm*B*PCOhigh^n)/(1+(B*PCOhigh^n));
qL = (qm*B*PCOlow^n)/(1+(B*PCOlow^n));
Mbed = (FCO)*300/((qH - qL)*fL); %mass of 1 PSA bed
VPSA =  Mbed/795/.8; %volume of 1 PSA vessel (m^3)
LPSA = (VPSA*64/pi)^(1/3);
DPSA = .25*LPSA;
ProductH2 = 0.9*FH2; %Hydrogen production rate (mol/s)
Exhaust = FCO2; %All CO2 goes into exhaust
FH2_out = 0.1*FH2;
FCO_out = FCO;
FCH4_out = FCH4;
Capacity = ProductH2/1033.399; %Fraction of desired H2 product
err = abs(Capacity-1);
Capacity2 = Capacity;
iter = iter + 1;
end
%% Energy Balances for Costing
Reactor1W = 1000*(-241.83*(FH2O_2 - FH2O_1) - 393.51*(FCO2_2 - FCO2_1)...
        - 110.53*(FCO_2 - FCO_1) - 74.87*(FCH4_2 - FCH4_1)); %Energy to heat reactor 1 (W)
Reactor2W = -1000*(-241.83*(FH2O_3 - FH2O_2) - 393.51*(FCO2_3 - FCO2_2)...
        - 110.53*(FCO_3 - FCO_2) - 74.87*(FCH4_3 - FCH4_2)); %Energy to cool reactor 1 (W)
PreHeaterW = (T1 - 298)*(35.69*FCH4_1 + 43.7934*FCO2_1 + 29.24*FH2_1); %Energy to heat methane feed to T1 (W)
BoilerW = 40700*FH2O_1 + 75.6*(212.5-80.7)*FH2O_1 + (380 - 212.5)*36.5*FH2O_1; %Energy to boil water in feed to steam at 380 Deg C (W)
ThermalW = BoilerW + PreHeaterW + Reactor1W;
FlashDrumW = (219-100)*(36.5*FH2O_3 + 35.69*FCH4_3 + 43.7934*FCO2_3 + 29.24*FH2_3 + 29.15*FCO_3) + 40700*FH2O_3; %Energy to cool stream and condense water (W)
CompressorBHP = (3.3*10^(-5)/.254)*P*14.5038*(FT_4*8.314*10^(-5)*373/P)*60*35.3147*((30/P)^.254 - 1)/.85; %BHP of compressor before PSA
CoolingWaterGJ = (FlashDrumW + Reactor2W)*3600*24*350/(10^9); %GJ of cooling water per year (GJ/year)
PumpW = 50150/1107.8*F2; %Pump duty (W)
%% Reactor Sizing Parameters
pcat = 1000; %Bulk density of catalyst (kg/m^3)
Vtot1 = W1/pcat; %Total Reactor 1 volume (m^3)
Vtot2 = W2/pcat;%Total Reactor 2 volume (m^3)
L1 = (Vtot1*20^2/pi)^(1/3); %Reactor 1 length (m) (assumes 10:1 aspect ratio)
Dr1 = .1*L1; %Reactor 1 Diameter (m) (assumes 10:1 aspect ratio)
L2 = (Vtot2*20^2/pi)^(1/3); %Reactor 2 length (m) (assumes 10:1 aspect ratio)
Dr2 = .1*L2; %Reactor 2 Diameter (m) (assumes 10:1 aspect ratio)
%% Variable Costs
H2Energy = 1000*286*FH2_out; %Energy released by burning excess H2 (W)
COEnergy = 1000*294*FCO_out; %Energy released by burning excess CO (W)
CH4Energy = 1000*890*FCH4_out; %Energy released by burning excess CH4 (W)
EffluentW = H2Energy + COEnergy + CH4Energy; %Total energy from burning all effluent (W)
EffluentCO2 = FCO_out + FCH4_out; %Total CO2 from burning all effluent (mol/s)
EnergyRatio = ThermalW/EffluentW; %Decides if more methane is needed for furnaces
if EnergyRatio < 1
    FurnaceCO2 = EnergyRatio*EffluentCO2;
    ThermalCost = 0; % Cost of thermal energy for furnaces ($/year)
else
    MethaneBurned = (ThermalW - EffluentW)/890000; %Total Methane burned in furnaces (mol/s) (discounts energy from burned effluent)
    FurnaceCO2 = MethaneBurned + EffluentCO2; %CO2 Produced by furnaces (mol/s)
    ThermalCost = (ThermalW - EffluentW)/1000000*10*24*350; % Cost of thermal energy for furnaces ($/year)
end
TotalElectricityW = CompressorBHP*745.7 + PumpW; %Total Electricity Usage (W)
ElectricityCO2 = 1.512*10^(-12)*TotalElectricityW; %CO2 Produced from Elec. (mol/s)
ElectricityCO2Cost = 8.8*10^(-3)*ElectricityCO2*3600*24*350; %Additional CO2 Cost from heaters ($/year)
FurnaceCO2Cost = 8.8*10^(-3)*FurnaceCO2*3600*24*350; %Additional CO2 Cost from heaters ($/year)
CO2Cost = 4.4*10^(-3)*Exhaust*3600*24*350 + FurnaceCO2Cost + ElectricityCO2Cost; %Total cost of sequestration ($/year)
FeedCost = .001879*F1*3600*24*350; %Cost of Methane Feed ($/year)
SteamFeedCost = .05*18.01/1000000*F2*3600*24*350; %Cost of Water Feed @ $0.05 per MT ($/year)
H2Revenue = H2Price/1000000*2.016*ProductH2*3600*24*350; %Income for H2 ($/year)
CoolingWaterCost = CoolingWaterGJ*2; %Cost of all cooling water ($/year)
CompressorElectricity = CompressorBHP*745.7/1000000*50*24*350; %Electricity powering compressor 1 ($/year)
PumpElectricity = PumpW/1000000*50*24*350; %Water pump electricity cost
WasteWaterCost = FH2O_out*3600*24*350*18.01/1000000*1.25; %Wastewater treatment cost ($/year)
AllCosts = [ThermalCost WasteWaterCost CO2Cost FeedCost SteamFeedCost CoolingWaterCost CompressorElectricity PumpElectricity]; %Vector of all variable costs
%% Fixed Costs
Fp = 1;
Fm = 1;
Fc = Fm*Fp;
PreHeaterBTU = PreHeaterW*3.412/1000000;
Reactor1BTU = Reactor1W*3.412/1000000;
BoilerBTU = BoilerW*3.412/1000000;
Reactor1IC = (1650/280)*(101.9*((Dr1*3.28)^1.066)*((L1*3.28)^.82))*(2.18+Fc);
Reactor2IC = (1650/280)*(101.9*((Dr2*3.28)^1.066)*((L2*3.28)^.82))*(2.18+Fc);
FeedFurnaceIC = (1650/280)*(5.52*10^3)*(PreHeaterBTU^0.85)*(1.27+Fc); % Pre-Reactor Furnace Installed Cost
BoilerIC = (1650/280)*(5.52*10^3)*(BoilerBTU^0.85)*(1.27+Fc);
Reactor1FurnaceIC = (1650/280)*(5.52*10^3)*(Reactor1BTU^0.85)*(1.27+Fc); % Reactor Furnace Installed Cost
CompressorIC = (1650/280)*517.5*CompressorBHP^.82*(2.11+Fc); %Installed cost of compressor 1
PSAIC = 4*(1650/280)*101.9*(DPSA*3.28)^1.066*(LPSA*3.28)^0.82*(3.18 + Fc) + 5*4*Mbed; %IC of 4 PSA tanks
HeatExchangerIC = (1650/280)*101.3*A1^.65*(2.29+Fc);
HeatExchangersIC = 16*HeatExchangerIC;
CatalystCost = (W1+W2)*30*3;
AllIC = [PSAIC Reactor1IC Reactor2IC FeedFurnaceIC BoilerIC Reactor1FurnaceIC CompressorIC HeatExchangersIC]; %Vector of all Installation Costs ($)
%% NPV Calculations
ISBL = sum(AllIC); %Total installation cost ($)
PR = H2Revenue; %AARON - annual product revenue ($/yr)
VC = sum(AllCosts);%AARON - annual cost of feedstock ($/yr)
FC = ISBL;
PBT = H2Revenue - sum(AllCosts);
 
Sigma = 6.8137;
a = 4.1107;
b = -0.3498;
c = 0.2936;
d = -2.06;
e = 3.4057;
f = 7.011;
g = 1.177;
h = 1.214;
k = 12;
 
TCI = h*FC;
ROIbt = PBT/TCI;
NPV = a*PBT+b*FC;
NPVP = c*ROIbt+d;
 
end
 
%% Functions for finding kinetic parameters
function k1 = k_1(T,R)
    k1 = 0.04*2.763*10^(-6)*exp((-240.1/R)*(1/T - 1/648));
end
function k2 = k_2(T,R)
    k2 = 0.04*.11337*exp((-67.13/R)*(1/T - 1/648));
end
function k3 = k_3(T,R)
    k3 = 0.04*3.289*10^(-7)*exp((-243.9/R)*(1/T - 1/648));
end
function k4 = k_4(T)
    k4 = 3.33*10^(-8)*exp(12.88-1855.5/T);
end
function K1 = K_1(T)
    K1 = exp(-26830/T + 30.114);
end
function K2 = K_2(T)
    K2 = exp(4400/T - 4.036);
end
function KCO = K_CO(T,R)
    KCO = 40.91*exp((70.65/R)*(1/T - 1/648));
end
function KH2 = K_H2(T,R)
    KH2 = 0.0296*exp((82.9/R)*(1/T - 1/648));
end
function KCH4 = K_CH4(T,R)
    KCH4 = 0.1791*exp((38.28/R)*(1/T - 1/823));
end
function KH2O = K_H2O(T,R)
    KH2O = 0.4152*exp((-88.68/R)*(1/T - 1/823));
end
function PSI = PSI_(KCO,KH2,KCH4,KH2O,PCO,PH2,PCH4,PH2O)
    PSI = 1 + KCO*PCO + KH2*PH2 + KCH4*PCH4 + (KH2O*PH2O/(PH2+.01));
end
function r1 = r_1(PSI,k1,K1,PH2,PCH4,PH2O,PCO)
    r1 = (k1/((PH2+.01)^2.5 * PSI^2))*(PCH4*PH2O - (PH2^3 * PCO/K1));
end
function r2 = r_2(PSI,k2,K2,PH2,PCO2,PH2O,PCO)
    r2 = (k2/((PH2+.01)*PSI^2))*(PCO*PH2O - (PH2*PCO2/K2));
end
function r3 = r_3(PSI,k3,K1,K2,PCH4,PH2,PCO2,PH2O)
    r3 = (k3/((PH2+.01)^3.5*PSI^2))*(PCH4*PH2O^2 - (PH2^4*PCO2/(K1*K2)));
end
function r4 = r_4(k4,K2,PCO,PH2O,PH2,PCO2)
    r4 = k4*PCO*PH2O*(1-(PH2*PCO2/(K2*PCO*PH2O)));
end