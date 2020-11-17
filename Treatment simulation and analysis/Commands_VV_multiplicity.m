%% Script to calculate the tumour burden for VV injections from 1 to 7


%% load fixed model parameters

% Viral therapy- The dosing strategy is loaded in the loop
%PA.ViralStartTime = 0;
PA.ViralOffset = 1; %Dosed at a maximum of once daily
PA.Viralkabs = 20; % absorption rate
PA.ViralAvailFrac = 1;

%Set baseline model parameters
PA.a1 =  1.662718177543053; %4T1 tumour growth parameters
PA.a2 = 1.444183287717442;%4T1 tumour growth parameters
PA.d1 =  0;
PA.d2 = 0.295260097294696;%4T1 tumour growth parameters
PA.d3 = PA.d2;
PA.tau = 33.7/24 - 1/PA.a2; % From Sato 2016 mean intermitotic time - the expected time in G1- gives a lower bound on a_2. 
PA.IntermitoticSD = 6.7/24;  % From Sato 2016 SD on intermitotic time

%Distribution Specific Parameters
PA.N = round(PA.tau.^2./PA.IntermitoticSD^2); % round((33.7./24)^2./(6.7./24));
PA.TransitRate = PA.N./PA.tau; % Transit rate across compartments
PA.d3Hat = PA.N./PA.tau.*(exp(PA.d3.*PA.tau./(PA.N+1))-1);
PA.d3HatR = PA.d3Hat;

%Viral parameters VV Strain
PA.kappaPrimer =  0.054058; %Infection rate [1/month]   ;
PA.deltaPrimer =  2.484864086006999; % Death rate of infected cells [1/month] ;
PA.alphaPrimer =  1.120771;  %Dying infected cell to virion conversion [virions/cell];
PA.omegaPrimer =  40.282073356365977;   %Virion decay rate [1/month].
PrimerDose = 1e2; %Amount of Primer (VV) administered

%Viral parameters VSV Strain (BOOSTER)
PA.kappaBooster =  0.065563768; %Infection rate [1/month]   ;
PA.deltaBooster =  11.003239918554174; % Death rate of infected cells [1/month] ;
PA.alphaBooster =  1.131130932;  %Dying infected cell to virion conversion [virions/cell];
PA.omegaBooster =  38.685206250628028;   %Virion decay rate [1/month].
BoosterDose = 0;%1e3; %Amount of Boosters (VSV) administered

%Strain independnent contact rate
PA.eta12 =  0.510538277701167;  % virion half effect contact rate, Strain independent

% Cytokine Parameters
PA.CprodHomeo = 0.00039863;   %Homeostatic cytokine production [ng/mL/month]
PA.CprodMax =  1.429574637713578;   % Maximal cytokine production rate [ng/mL/month]
PA.C12 = 0.739376299393775; % Half effect in cytokine production [cells/month]
PA.kel = 0.16139;   % elimination rate of cytokine [1/month]

%Immune Parameters
PA.kp = 9.23124604834137;    % contact rate with phagocytes [1/month]
PA.kq = 0.06483950327214;  % factor in denominator of Q phagocytosis term [Unitless] 
PA.ks = PA.kq; % factor in denominator of S phagocytosis term [Unitless]
PA.P12 = 0.000114983672183;    % Half effect in cytokine driven phagocyte production [ng/mL]
PA.gammaP = 0.35; % From Barrish 2017 PNAS elimination rate of phagocyte [1/month]
PA.Kcp = 4.6754;  % From Barrish 2017 PNAS conversion of cytokine into phagocyte [cells/ng/mL]
PA.KcpBooster = 3.081693813811469; % Kcp when VV has primed the tumour microenvironment and VSV is injected

%Resistant Parameters
PA.nu = 1e-10;  %  Mutation probability
PA.a1R = PA.a1;
PA.a2R = PA.a2;
PA.d1R = PA.d1;
PA.d2R = PA.d2;
PA.d3R = PA.d3;


%Immunotherapy (Turned off)
PA.AdminNumber = 1; %Need this for the code to run, even if no immunotherapy is actually administered
PA.StartTime = 0;
PA.Offset = 1; %Dosed every day 
PA.Admin = 0; % Administer no immunotherapy 
PA.Vol = 7; %volume of absorption
PA.kabs = 6.6311; % absorption rate
PA.AvailFrac =0.85; 

%% Load the patient-specific matrix of parameters
%Load the patient specific parameters
VirtualPatientParametersInput =load('Virtual_patient_parameters');
load('Initial_Cond'); %Load the 200 virtual patient initial conditions
B = struct2array(VirtualPatientParametersInput);
M = length(B(1,:)); %The number of patients

NumberOfPatients = 200;
TumourBurdenMatrix = zeros(NumberOfPatients,7); % matrix to store the tumour burden for each patient with jj priming doses given

dcount = 1;
%% Calculate tumour burden for each patient and number of primer doses

for ii = 1:200 % simulation for all patients
    %   Update the parameters for each virtual patient
    Param = B([1,2,3,4,6,21,22,23,26],ii+1); %load the varied parameters from the virtual patient and update PA.
     C = {Param(1), ... 
         Param(2), ...
         Param(3), ...
         Param(4), ...
         Param(5), ...
         Param(6), ...
         Param(7), ...
         Param(8), ...
         Param(9)};
    ParameterNames = {'a1' 'a2' 'd1' 'd2' 'tau' 'kp' 'kq' 'ks' 'Kcp'}; %The name of parameters to be varied
    %  Update model parameters
   [PA.(ParameterNames{1}),PA.(ParameterNames{2}),PA.(ParameterNames{3}),PA.(ParameterNames{4})...
        PA.(ParameterNames{5}),PA.(ParameterNames{6}),PA.(ParameterNames{7}),PA.(ParameterNames{8}),...
        PA.(ParameterNames{9})] = C{:}; %update the parameters for this run 
 
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
    
    for jj = 1:7 %possible VV doses
        PA.ViralAdminNumberStrainPrimer =jj; % maximal number of VV
        PA.ViralAdminNumberStrainBooster = 0;%1; % maximal number of VSV
        
        ViralTherapyDosesStrainPrimer = ones(1,jj); % Create a dose vector with jj VV dosess administered on consecutive days
        ViralTherapyDosesStrainBooster = ones(1,1); % Create a dose vector with jj VSV dosess administered on consecutive days
        
        PA.ViralStartTimePrimer = 0; %first VV injection daY
        PA.ViralStartTimeBooster = jj-1+7; %VSV administration day
        
        PA.ViralAdminStrainPrimer =  ViralTherapyDosesStrainPrimer.*PrimerDose; %Amount of VV administered
 
        PA.ViralAdminStrainBooster = ViralTherapyDosesStrainBooster.*BoosterDose; %Amount of VSV administered
        
        % Time interval
        tf = (50); % End time is 15 days after last priming dose.
        tstep = floor(tf./20);
        totaltime = [0 tf];
        
        [solTreat PAN] =  model_simulator(PA,totaltime,Initial_cond_day_6(ii,:));
       
        %calculating tumour burden from last primer
        TimeSeries = linspace(jj-1,jj-1+7+15,1001); %Create 1000 evenly space points between 0 and tf
        IntStep = TimeSeries(2)-TimeSeries(1); %Calculate the time step between points
        EvalSol = deval(solTreat,TimeSeries); %Evaluate the solution at the collocation points
        TrapFun = EvalSol(1,:)+ EvalSol(2,:)+ EvalSol(3,:)+ EvalSol(PAN+7,:)+ EvalSol(PAN+8,:)+ EvalSol(PAN+9,:)+ EvalSol(PAN+PAN+10,:)+ EvalSol(PAN+PAN+12,:); % Find the tumour burden at the collocation points.
        Tumour_dynam{ii}(jj,:) = TrapFun;
        
        % Calculate tumour burden: Area under the tumour curve between tf/2  and tf.
        TimeSeries = linspace(0,jj-1+7+15,1001); %Create 1000 evenly space points between 0 and tf
        IntStep = TimeSeries(2)-TimeSeries(1); %Calculate the time step between points
        EvalSol = deval(solTreat,TimeSeries); %Evaluate the solution at the collocation points
        TrapFun = EvalSol(1,:)+ EvalSol(2,:)+ EvalSol(3,:)+ EvalSol(PAN+7,:)+ EvalSol(PAN+8,:)+ EvalSol(PAN+9,:)+ EvalSol(PAN+PAN+10,:)+ EvalSol(PAN+PAN+12,:); % Find the tumour burden at the collocation points.
        TumourBurden_Final{dcount}(ii,jj) = TrapFun(end); 
        
    end
    PatientNumber = ii
end

%%

figure
boxplot([TumourBurden_Final{1}])
xlabel('Number of primers (D_p)')
ylabel('Tumour Burden')
set(gca,'FontSize',13)

