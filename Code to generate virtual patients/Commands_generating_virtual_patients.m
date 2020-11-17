%% This code generates virtual patients by sampling parameters from normal 
%distributions and checking their tumour growths lie within reasonable 
%standard deviation from the data

close all
clear all
format long

%% Input data (from Le Boeuf et al. 2010 "Syndergistic Interaction Between On
LeBoeufSamplingTime = [6,8,10,13,15,18]; 
LeBoeufImmuneData = [0,4,13,18,37,55];
ErrorThreshold = 0.5; % Acceptable amount of error: ResidualGenerated/ResidualFitted \in ( 1-Threshold,1+Threshold)
tf = LeBoeufSamplingTime(end);
totaltime = [0 tf];

%% Set baseline model parameters
% Growth parameters
PA.a1 =  1.662718177543053; %Transit Rate Quiescent [1/month]
PA.a2 = 1.444183287717442; %Transit Rate G1 [1/month]
PA.d1 =  0; %Death Rate Quiescent [1/month]
PA.d2 =  0.295260097294696; %Death Rate G1 [1/month]
PA.d3 = PA.d2; %Death rate Mitotic [1/month]
%Delay parameters
PA.tau = 33.7/24-1/PA.a2; % From Sato 2016 mean intermitotic time - the expected time in G1- gives a lower bound on a_2. 
PA.IntermitoticSD = 6.7/24; % From Sato 2016 SD on intermitotic time

%Distribution Specific Parameters
PA.N = round(PA.tau.^2./PA.IntermitoticSD^2); 
PA.TransitRate = PA.N./PA.tau; % Transit rate across compartments
PA.d3Hat = PA.N./PA.tau.*(exp(PA.d3.*PA.tau./(PA.N+1))-1);
PA.d3HatR = PA.d3Hat;

%Viral parameters
PA.kappa =  3.534412642851458; %Infection rate [1/month]   ;
PA.delta =  4.962123414821151; % Death rate of infected cells [1/month] ;
PA.alpha =  0.008289097649957;  %Dying infected cell to virion conversion [virions/cell];
PA.omega =  9.686308020782763;   %Virion decay rate [1/month].
PA.eta12 =  0.510538277701167;  % virion half effect contact rate

% Cytokine Parameters
PA.CprodHomeo = 0.00039863;   %Homeostatic cytokine production [ng/mL/month]
PA.CprodMax =  1.429574637713578;   % Maximal cytokine production rate [ng/mL/month]
PA.C12 = 0.739376299393775; % Half effect in cytokine production [cells/month]
PA.kel = 0.16139;   % elimination rate of cytokine [1/month]

%Immune Parameters
PA.kp =  9.23124604834137; %Infection rate [1/month]   ;
PA.kq =  0.064839; % Death rate of infected cells [1/month] ;
PA.ks = PA.kq; % factor in denominator of S phagocytosis term [Unitless]
PA.P12 = 0.0001149;    % Half effect in cytokine driven phagocyte production [ng/mL]
PA.gammaP = 0.35; % From Barrish 2017 PNAS elimination rate of phagocyte [1/month]
PA.Kcp = 4.6754; % From Barrish 2017 PNAS conversion of cytokine into phagocyte [cells/ng/mL]

% Immune Steady State
PA.CStar = PA.CprodHomeo./PA.kel;
PA.PStar = (1./PA.gammaP).*(PA.Kcp.*PA.CStar./(PA.P12+PA.CStar));

%Resistant Parameters
PA.nu = 1e-10; % Mutation percentage
PA.a1R = PA.a1;
PA.a2R = PA.a2;
PA.d1R = PA.d1;
PA.d2R = PA.d2;
PA.d3R = PA.d3;
PA.kappaR = PA.kappa;
 
%% Set up sampling procedure
PA.TotalCells = 1e5; 
PA.DeathTime = 2 ;

SizePA = length(struct2array(PA)); %How many parameters are there

%Parameter set up
ParameterNames = {'a1' 'a2' 'd1' 'd2' 'tau' 'kp' 'kq' 'ks' 'Kcp'}; %The name of parameters to be varied
ParameterHomeostasis = [PA.a1 ;PA.a2 ; PA.d1; PA.d2; PA.tau; PA.kp; PA.kq; PA.ks; PA.Kcp];
N = 0; % Number of viable patients discovered.
n = 0; %Number of non-viable patients
M = 200; %Number of desired virtual patients

VirtualPatientParameters = zeros(SizePA,M); %array to save the parameters for each virtual patient
ResidualGeneratedErrorVector = zeros(1,M+1); % Vector to record the Generated Error

A = repmat(ParameterHomeostasis,[1,100*(M+1)]);
S =  repmat((0.04./2.5).*ParameterHomeostasis,[1,100*(M+1)]);
B = normrnd(A,S,9,100*(M+1)); % Matrix of 100xM normally distributed parameters with mean homeostatic parameter and standard deviation S.

%% Solve the fitted system to calcuate the Fitted Residual Error
%Calculate the cell cycle duration
TotalTime = 1/PA.a1 +1/(PA.a2+PA.d2)+PA.tau;
% Initial Conditions
QIC = (1/PA.a1./TotalTime).*PA.TotalCells.*(1-PA.nu); %100;
SIC = (1/(PA.a2+PA.d2)./TotalTime).*PA.TotalCells.*(1-PA.nu);  %100;
TCIC = (PA.tau./TotalTime).*PA.TotalCells.*(1-PA.nu).*ones(1,PA.N)./PA.N; %Transit compartment ICs
NCIC = (PA.tau./TotalTime).*PA.TotalCells.*(1-PA.nu);
IIC = 0;
VIC = 0;
CIC =  PA.CprodHomeo./PA.kel;
PIC =   PA.Kcp.*CIC./((PA.P12+CIC).*PA.gammaP);
% Resistant Strain ICs
RIC =   (1/PA.a1./TotalTime).*PA.TotalCells.*(PA.nu); %
RSIC =   (1/(PA.a2+PA.d2)./TotalTime).*PA.TotalCells.*(PA.nu); %
ResistantTCIC =   (PA.tau./TotalTime).*PA.TotalCells.*(PA.nu).*ones(1,PA.N)./PA.N; 
ResistantTotalCellsIC =   (PA.tau./TotalTime).*PA.TotalCells.*(PA.nu);

InitialConditionsFit = [QIC,SIC,IIC,VIC,TCIC,CIC,PIC,NCIC,RIC,RSIC,ResistantTCIC,ResistantTotalCellsIC];

[solFitted] =  WaresDistributedImmunityResistantSolver(totaltime,InitialConditionsFit,PA);
LeBoeufFittedDataErrorVec = deval(solFitted,LeBoeufSamplingTime,1)+ deval(solFitted,LeBoeufSamplingTime,2)+ deval(solFitted,LeBoeufSamplingTime,PA.N+7)...
                             + deval(solFitted,LeBoeufSamplingTime,PA.N+8)+ deval(solFitted,LeBoeufSamplingTime,PA.N+9)+ deval(solFitted,LeBoeufSamplingTime,(PA.N+PA.N+10)) ;
LeBoeufFittedDataRelativeErrorVec = 100.*(LeBoeufFittedDataErrorVec(:)-LeBoeufFittedDataErrorVec(1))./LeBoeufFittedDataErrorVec(1);
FittedResidualError = sqrt(dot(LeBoeufFittedDataRelativeErrorVec-LeBoeufImmuneData',LeBoeufFittedDataRelativeErrorVec-LeBoeufImmuneData'));

Bad_params = [];

while N <   M+1 
    %Sample from Parameter space
    Param = B(:,n+N+1);
    C = {Param(1), ... 
      Param(2), ...
      Param(3), ...
      Param(4), ...
      Param(5), ...
      Param(6), ...
      Param(7), ...
      Param(8), ...
      Param(9) };

    %Update model parameters
    [PA.(ParameterNames{1}),PA.(ParameterNames{2}),PA.(ParameterNames{3}),PA.(ParameterNames{4})...
     PA.(ParameterNames{5}),PA.(ParameterNames{6}),PA.(ParameterNames{7}),PA.(ParameterNames{8}),PA.(ParameterNames{9})] = C{:}; %update the parameters for this run 
 
    % Resistant parameters
    PA.a1R = PA.a1;
    PA.a2R = PA.a2;
    PA.d1R = PA.d1;
    PA.d2R = PA.d2;
    PA.d3R = PA.d3;
    
    PA.TransitRate = PA.N./PA.tau; %must recalculate transit rate, as the delay varies
    PA.CStar = PA.CprodHomeo./PA.kel;
    PA.PStar = (1./PA.gammaP).*(PA.Kcp.*PA.CStar./(PA.P12+PA.CStar));
    PA.d3Hat = PA.N./PA.tau.*(exp(PA.d3.*PA.tau./(PA.N+1))-1);
    PA.d3HatR = PA.d3Hat;
    
    %% Calculate the initial conditions
    %Calculate the cell cycle duration
    TotalTime = 1/PA.a1 +1/(PA.a2+PA.d2)+PA.tau;
    % Initial Conditions
    QIC = (1/PA.a1./TotalTime).*PA.TotalCells.*(1-PA.nu); %100;
    SIC = (1/(PA.a2+PA.d2)./TotalTime).*PA.TotalCells.*(1-PA.nu);  %100;
    TCIC = (PA.tau./TotalTime).*PA.TotalCells.*(1-PA.nu).*ones(1,PA.N)./PA.N; %Transit compartment ICs
    NCIC = (PA.tau./TotalTime).*PA.TotalCells.*(1-PA.nu);
    IIC = 0;
    VIC = 0;
    CIC =  PA.CprodHomeo./PA.kel;
    PIC =  PA.Kcp.*CIC./((PA.P12+CIC).*PA.gammaP);
    % Resistant Strain ICs
    RIC =   (1/PA.a1./TotalTime).*PA.TotalCells.*(PA.nu); %
    RSIC =   (1/(PA.a2+PA.d2)./TotalTime).*PA.TotalCells.*(PA.nu); %
    ResistantTCIC =   (PA.tau./TotalTime).*PA.TotalCells.*(PA.nu).*ones(1,PA.N)./PA.N; 
    ResistantTotalCellsIC =   (PA.tau./TotalTime).*PA.TotalCells.*(PA.nu);

    InitialConditionsTest = [QIC,SIC,IIC,VIC,TCIC,CIC,PIC,NCIC,RIC,RSIC,ResistantTCIC,ResistantTotalCellsIC];

    SEM = [1,3,4,6,5];
    threeSTDfromfit = SEM.*sqrt(10)*2;

    %%  Solve the Distributed DDE system

    [solTest] =  WaresDistributedImmunityResistantSolver(totaltime,InitialConditionsTest,PA);
    LeBoeufDataErrorVecMouse =  deval(solTest,LeBoeufSamplingTime,1)+ deval(solTest,LeBoeufSamplingTime,2)+ deval(solTest,LeBoeufSamplingTime,PA.N+7)...
                                 + deval(solTest,LeBoeufSamplingTime,PA.N+8)+ deval(solTest,LeBoeufSamplingTime,PA.N+9)+ deval(solTest,LeBoeufSamplingTime,(PA.N+PA.N+10)) ;
    LeBoeufMouseDataRelativeErrorVec = 100.*(LeBoeufDataErrorVecMouse(:)-LeBoeufDataErrorVecMouse(1))./LeBoeufDataErrorVecMouse(1);
    ResidualGeneratedError = sqrt(dot(LeBoeufMouseDataRelativeErrorVec-LeBoeufImmuneData',LeBoeufMouseDataRelativeErrorVec-LeBoeufImmuneData'));
    ResidualErrorComparison = ResidualGeneratedError/FittedResidualError;


    RE = abs(LeBoeufMouseDataRelativeErrorVec(2:end)-LeBoeufImmuneData(2:end)');
    errorvec = RE' - threeSTDfromfit;

    gradient = LeBoeufDataErrorVecMouse(end)-LeBoeufDataErrorVecMouse(end-1);
    if 1-ErrorThreshold < ResidualErrorComparison && isempty(find(errorvec>0))==1 %&& gradient>0 %&& LeBoeufMouseDataRelativeErrorVec(end)>50%ResidualErrorComparison < 2.25  %1+ErrorThreshold %Test to see if ResidualGeneratedError/FittedResidualError is acceptable
         N = N+1
         ResidualGeneratedErrorVector(N+1) = ResidualErrorComparison;
         VirtualPatientParameters(:,N+1) = struct2array(PA)'; %Save the virtual patient parameters
    else
        Bad_params = [Bad_params; [ResidualErrorComparison Param([1,2,4:9])']];
        n = n+1
       % ResidualErrorComparison
    end

end
TotalNumberofRuns = N+n

save('Virtual_patient_parameters.mat','VirtualPatientParameters')

%% Create plots
for ii = 1:M+1
  	Param = VirtualPatientParameters([1,2,3,4,6,21,22,23,26],ii+1); %load the varied parameters from the virtual patient and update PA.
    CTest = {Param(1), ... 
          Param(2), ...
          Param(3), ...
          Param(4), ...
          Param(5), ...
          Param(6), ...
          Param(7), ...
          Param(8), ...
          Param(9) };

    %Update model parameters
    [PA.(ParameterNames{1}),PA.(ParameterNames{2}),PA.(ParameterNames{3}),PA.(ParameterNames{4})...
     PA.(ParameterNames{5}),PA.(ParameterNames{6}),PA.(ParameterNames{7}),PA.(ParameterNames{8}),PA.(ParameterNames{9})] = CTest{:}; %update the parameters for this run 
 
 
    % Resistant parameters
    PA.a1R = PA.a1;
    PA.a2R = PA.a2;
    PA.d1R = PA.d1;
    PA.d2R = PA.d2;
    PA.d3R = PA.d3;
    
    PA.TransitRate = PA.N./PA.tau; %must recalculate transit rate, as the delay varies
    PA.CStar = PA.CprodHomeo./PA.kel;
    PA.PStar = (1./PA.gammaP).*(PA.Kcp.*PA.CStar./(PA.P12+PA.CStar));
    PA.d3Hat = PA.N./PA.tau.*(exp(PA.d3.*PA.tau./(PA.N+1))-1);
    PA.d3HatR = PA.d3Hat;
    
    %% Calculate the initial conditions
    %Calculate the cell cycle duration
    TotalTime = 1/PA.a1 +1/(PA.a2+PA.d2)+PA.tau;
    % Initial Conditions
    QIC = (1/PA.a1./TotalTime).*PA.TotalCells.*(1-PA.nu); %100;
    SIC = (1/(PA.a2+PA.d2)./TotalTime).*PA.TotalCells.*(1-PA.nu);  %100;
    TCIC = (PA.tau./TotalTime).*PA.TotalCells.*(1-PA.nu).*ones(1,PA.N)./PA.N; %Transit compartment ICs
    NCIC = (PA.tau./TotalTime).*PA.TotalCells.*(1-PA.nu);
    IIC = 0;
    VIC = 0;
    CIC =  PA.CprodHomeo./PA.kel;
    PIC =   PA.Kcp.*CIC./((PA.P12+CIC).*PA.gammaP);
    % Resistant Strain ICs
    RIC =   (1/PA.a1./TotalTime).*PA.TotalCells.*(PA.nu); %
    RSIC =   (1/(PA.a2+PA.d2)./TotalTime).*PA.TotalCells.*(PA.nu); %
    ResistantTCIC =   (PA.tau./TotalTime).*PA.TotalCells.*(PA.nu).*ones(1,PA.N)./PA.N; 
    ResistantTotalCellsIC =   (PA.tau./TotalTime).*PA.TotalCells.*(PA.nu);

    InitialConditionsCheck = [QIC,SIC,IIC,VIC,TCIC,CIC,PIC,NCIC,RIC,RSIC,ResistantTCIC,ResistantTotalCellsIC];

    [solPlot] =  WaresDistributedImmunityResistantSolverTest(totaltime,InitialConditionsCheck,PA); 

    InitializeValueMousePlot = deval(solPlot,LeBoeufSamplingTime,1)+ deval(solPlot,LeBoeufSamplingTime,2)+ deval(solPlot,LeBoeufSamplingTime,PA.N+7)...
                                 + deval(solPlot,LeBoeufSamplingTime,PA.N+8)+ deval(solPlot,LeBoeufSamplingTime,PA.N+9)+ deval(solPlot,LeBoeufSamplingTime,PA.N+PA.N+10); 

   Fig1 = figure(1);
   hold on
   g1 = plot(solPlot.x(:),100.*( solPlot.y(1,:)+solPlot.y(2,:)+solPlot.y(PA.N+7,:)+solPlot.y(PA.N+8,:)+solPlot.y(PA.N+9,:)+solPlot.y(PA.N+PA.N+10,:)- InitializeValueMousePlot(1))./(InitializeValueMousePlot(1)),'LineWidth',2,'Color','b','LineStyle','-');
   xlabel('Time (days)','FontSize',15)
   ylabel('Relative Tumour Size','FontSize',25,'Interpreter','latex','FontSize',15)
end
figure(1)
InitializeValue = LeBoeufFittedDataErrorVec;
g1 = plot(solFitted.x(:),100.*( solFitted.y(1,:)+solFitted.y(2,:)+solFitted.y(PA.N+7,:)+solFitted.y(PA.N+8,:)+solFitted.y(PA.N+9,:)+solFitted.y(PA.N+PA.N+10,:)- InitializeValue(1))./(InitializeValue(1)),'LineWidth',2,'Color','k','LineStyle','-');
hold on
scatter(LeBoeufSamplingTime,LeBoeufImmuneData,40,'r','*');

%% Parameter Test Solver
function [sol] = WaresDistributedImmunityResistantSolver(totaltime,IC,PA) %DDE model without therapy
% opts = ddeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',1e-2,'Events',@EventsViralMonths1);
opts = odeset('RelTol',1e-8,'AbsTol',1e-8,'MaxStep',1e-2,'Events',@EventsViralMonths1);
sol = ode15s(@ViralOncologyParameterFit,totaltime,IC,opts);
%opts = ddeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',1e-2,'Events',@EventsViralMonths1);
%sol = ddesd_f5(@Wares2012DistributedImmunity,@(t,y) DelayWares1(t,y,PA) ,IC,totaltime,opts);
function dydt = ViralOncologyParameterFit(t,y,Z);
%Quiescent cells
dydt(1) = 2.*(1-PA.nu).*PA.TransitRate.*y(PA.N+4)-(PA.a1+PA.d1+ psiQ(y(PA.N+6),y(1),PA) ).*y(1); 
%G1 Cells
dydt(2) = PA.a1.*y(1)-( PA.a2+PA.d2+ PA.kappa.*Infection(y(4),PA)+ psiS(y(PA.N+6),y(2),PA) ).*y(2); 
%Infected Cells
dydt(3) =  0; % -PA.delta.*y(3)+ PA.kappa.*Infection(y(4),PA).*(y(2)+y(PA.N+7)+y(PA.N+9)+y(PA.N+PA.N+10) ); % PA.kappa.*Infection(y(1),y(2)+Sbar,y(3),y(4),PA).*(y(2)+Sbar); % Infected cells
%Virions
dydt(4) =   0; % PA.alpha.*PA.delta.*y(3)-PA.omega.*y(4)- PA.kappa.*Infection(y(4),PA).*(y(2)+y(PA.N+7)+y(PA.N+9)+y(PA.N+PA.N+10)); % PA.kappa.*Infection(y(1),y(2)+Sbar,y(3),y(4),PA ).*(y(2)+Sbar); % Virions
%Writing ODE for first transit compartment
 dydt(5) = PA.a2.*y(2) - PA.TransitRate.*y(5) - ( PA.d3Hat+PA.kappa.*Infection(y(4),PA) + psiS(y(PA.N+6),y(5),PA) ).*y(5);
for jj = 6:PA.N+4
    dydt(jj) = PA.TransitRate.*(y(jj-1)-y(jj)) - ( PA.d3Hat +PA.kappa.*Infection(y(4),PA) + psiS(y(PA.N+6),y(jj),PA) ).*y(jj); %Transit compartment ODEs
end
%Immune cytokine 
dydt(PA.N+5) =   CProd( psiS(y(PA.N+6),y(2),PA).*(y(2)+y(PA.N+7))+PA.delta.*y(3)+psiQ(y(PA.N+6),y(1),PA).*y(1),PA) - PA.kel.*y(PA.N+5); 
%Phagocytes
dydt(PA.N+6) =   PA.Kcp.*y(PA.N+5)./(PA.P12+y(PA.N+5)) - PA.gammaP.*y(PA.N+6); 
%ODE for total number of cells in cell cycle
dydt(PA.N+7) = PA.a2.*y(2) - (PA.d3Hat+PA.kappa.*Infection(y(4),PA)+psiS(y(PA.N+6),y(2),PA) ).*y(PA.N+7) - (PA.TransitRate./PA.a2).*y(PA.N+4);
% Resistant compartments
%Resistant Quiescent
dydt(PA.N+8) =   2.*PA.nu.*PA.TransitRate.*y(PA.N+4) + 2.*PA.TransitRate.*y(PA.N+PA.N+9)-(PA.a1R+PA.d1R).*y(PA.N+8); %Resistant quiescence DDE
%Resistant G1
dydt(PA.N+9) =   PA.a1R.*y(PA.N+8)-(PA.a2R+PA.d2R+ PA.kappa.*Infection(y(4),PA) ).*y(PA.N+9); %Susceptible resistant cells
%Resistant First transit
dydt(PA.N+10) =  PA.a2R.*y(PA.N+9) - PA.TransitRate.*y(PA.N+10) - (PA.d3HatR+PA.kappa.*Infection(y(4),PA)).*y(PA.N+10); %Susceptible resistant first compartment
for jj = PA.N+11:PA.N+PA.N+9
    dydt(jj) =  PA.TransitRate.*(y(jj-1)-y(jj)) - (PA.d3HatR +PA.kappa.*Infection(y(4),PA)  ).*y(jj); %Resistant Transit compartment ODEs
end
%DE for total resistant cells
dydt(PA.N+PA.N+10) =   PA.a2.*y(PA.N+9) - (PA.d3Hat+PA.kappa.*Infection(y(4),PA) ).*y(PA.N+PA.N+10) - (PA.TransitRate./PA.a2).*y(PA.N+PA.N+9);
dydt = dydt';
end
function [value,isterminal,direction] = EventsViralMonths1(t,y,Z)

value(1) = y(1)+ y(2) +y(PA.N+7)  - 2.*PA.TotalCells;  %What we are setting to 0, the tumour doubling size (different initial conditions, but close enough

isterminal(1) = 0;   % 1 = End the integration
direction(1) = 0;   % Positive direction only

value(2)= y(1)+ y(2)+ y(PA.N+7) - PA.TotalCells.*2^(13+PA.DeathTime);  %What we are setting to 0 
isterminal(2) = 1;   % 1 = End the integration
direction(2) = 0;   % Positive direction only

end
end

function [sol] = WaresDistributedImmunityResistantSolverTest(totaltime,IC,PA) %DDE model without therapy
% opts = ddeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',1e-2,'Events',@EventsViralMonths1);
opts = odeset('RelTol',1e-8,'AbsTol',1e-8,'MaxStep',1e-2,'Events',@EventsViralMonths1);
sol = ode15s(@ViralOncologyParameterFit,totaltime,IC,opts);
%opts = ddeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',1e-2,'Events',@EventsViralMonths1);
%sol = ddesd_f5(@Wares2012DistributedImmunity,@(t,y) DelayWares1(t,y,PA) ,IC,totaltime,opts);
function dydt = ViralOncologyParameterFit(t,y,Z);
%Quiescent cells
dydt(1) = 2.*(1-PA.nu).*PA.TransitRate.*y(PA.N+4)-(PA.a1+PA.d1+ psiQ(y(PA.N+6),y(1),PA) ).*y(1); 
%G1 Cells
dydt(2) = PA.a1.*y(1)-( PA.a2+PA.d2+ PA.kappa.*Infection(y(4),PA)+ psiS(y(PA.N+6),y(2),PA) ).*y(2); 
%Infected Cells
dydt(3) =  0; % -PA.delta.*y(3)+ PA.kappa.*Infection(y(4),PA).*(y(2)+y(PA.N+7)+y(PA.N+9)+y(PA.N+PA.N+10) ); % PA.kappa.*Infection(y(1),y(2)+Sbar,y(3),y(4),PA).*(y(2)+Sbar); % Infected cells
%Virions
dydt(4) =   0; % PA.alpha.*PA.delta.*y(3)-PA.omega.*y(4)- PA.kappa.*Infection(y(4),PA).*(y(2)+y(PA.N+7)+y(PA.N+9)+y(PA.N+PA.N+10)); % PA.kappa.*Infection(y(1),y(2)+Sbar,y(3),y(4),PA ).*(y(2)+Sbar); % Virions
%Writing ODE for first transit compartment
 dydt(5) = PA.a2.*y(2) - PA.TransitRate.*y(5) - ( PA.d3Hat+PA.kappa.*Infection(y(4),PA) + psiS(y(PA.N+6),y(5),PA) ).*y(5);
for jj = 6:PA.N+4
    dydt(jj) = PA.TransitRate.*(y(jj-1)-y(jj)) - ( PA.d3Hat +PA.kappa.*Infection(y(4),PA) + psiS(y(PA.N+6),y(jj),PA) ).*y(jj); %Transit compartment ODEs
end
%Immune cytokine 
dydt(PA.N+5) =   CProd( psiS(y(PA.N+6),y(2),PA).*(y(2)+y(PA.N+7))+PA.delta.*y(3)+psiQ(y(PA.N+6),y(1),PA).*y(1),PA) - PA.kel.*y(PA.N+5); 
%Phagocytes
dydt(PA.N+6) =   PA.Kcp.*y(PA.N+5)./(PA.P12+y(PA.N+5)) - PA.gammaP.*y(PA.N+6); 
%ODE for total number of cells in cell cycle
dydt(PA.N+7) = PA.a2.*y(2) - (PA.d3Hat+PA.kappa.*Infection(y(4),PA)+psiS(y(PA.N+6),y(2),PA) ).*y(PA.N+7) - (PA.TransitRate./PA.a2).*y(PA.N+4);
% Resistant compartments
%Resistant Quiescent
dydt(PA.N+8) =   2.*PA.nu.*PA.TransitRate.*y(PA.N+4) + 2.*PA.TransitRate.*y(PA.N+PA.N+9)-(PA.a1R+PA.d1R).*y(PA.N+8); %Resistant quiescence DDE
%Resistant G1
dydt(PA.N+9) =   PA.a1R.*y(PA.N+8)-(PA.a2R+PA.d2R+ PA.kappa.*Infection(y(4),PA) ).*y(PA.N+9); %Susceptible resistant cells
%Resistant First transit
dydt(PA.N+10) =  PA.a2R.*y(PA.N+9) - PA.TransitRate.*y(PA.N+10) - (PA.d3HatR+PA.kappa.*Infection(y(4),PA)).*y(PA.N+10); %Susceptible resistant first compartment
for jj = PA.N+11:PA.N+PA.N+9
    dydt(jj) =  PA.TransitRate.*(y(jj-1)-y(jj)) - (PA.d3HatR +PA.kappa.*Infection(y(4),PA)  ).*y(jj); %Resistant Transit compartment ODEs
end
%DE for total resistant cells
dydt(PA.N+PA.N+10) =   PA.a2.*y(PA.N+9) - (PA.d3Hat+PA.kappa.*Infection(y(4),PA) ).*y(PA.N+PA.N+10) - (PA.TransitRate./PA.a2).*y(PA.N+PA.N+9);
dydt = dydt';
end
function [value,isterminal,direction] = EventsViralMonths1(t,y,Z)

value(1) = y(1)+ y(2) +y(PA.N+7)  - 2.*PA.TotalCells;  %What we are setting to 0, the tumour doubling size (different initial conditions, but close enough

isterminal(1) = 0;   % 1 = End the integration
direction(1) = 0;   % Positive direction only

value(2)= y(1)+ y(2)+ y(PA.N+7) - PA.TotalCells.*2^(13+PA.DeathTime);  %What we are setting to 0 
isterminal(2) = 1;   % 1 = End the integration
direction(2) = 0;   % Positive direction only

end
end


function g = psiQ(P,Q,PA)
% clearance of quiescent cells by immune system
 g = P.*(PA.kp./(1+PA.kq.*Q));
% g = PA.kp.*P;
end

function h = psiS(P,S,PA)
% clearance of susceptible cells by immune system
h =  P.*(PA.kp./(1+PA.ks.*S));
% h = P.*PA.kp;
end

function y = CProd(a,PA)
%Cytokine Production
y = PA.CprodHomeo+ (PA.CprodMax-PA.CprodHomeo).*(a./(PA.C12+a));
end

function f = Infection(V,PA) %Contact function
    if V > 1e-10
        f = V./(PA.eta12+V);
    else
        f = 0;
    end
end

% Immunotherapy dosing function
function DoseWares = Dose(PA,t);
DoseVec = zeros(1,PA.AdminNumber);
TAdmin = PA.StartTime+(0:PA.AdminNumber-1).*PA.Offset-t;
for nn = 1:PA.AdminNumber
    if TAdmin(nn) < 0
    DoseVec(nn) = (PA.AvailFrac.*PA.kabs.*PA.Admin(nn))./PA.Vol.*exp(PA.kabs.*TAdmin(nn));
    else
        DoseVec(nn) = 0;
    end
end
DoseWares = ones(1,PA.AdminNumber)*DoseVec';
end

% Viral dosing function
function DoseViral = ViralDose(PA,t);
ViralDoseVec = zeros(1,PA.ViralAdminNumber);
TAdminViral = PA.ViralStartTime+(0:PA.ViralAdminNumber-1).*PA.ViralOffset-t;
for nn = 1:PA.ViralAdminNumber
    if TAdminViral(nn) < 0
    ViralDoseVec(nn) = (PA.Viralkabs.*PA.ViralAdmin(nn).*PA.ViralAvailFrac)./PA.Vol.*exp(PA.Viralkabs.*TAdminViral(nn));
    else
        ViralDoseVec(nn) = 0;
    end
end
DoseViral = ones(1,PA.ViralAdminNumber)*ViralDoseVec';
end

function d = DelayWares1(t,y,PA)
%This function sets up the delay vectors necessary for the DDE solver.
d = [t-PA.tau];            
end