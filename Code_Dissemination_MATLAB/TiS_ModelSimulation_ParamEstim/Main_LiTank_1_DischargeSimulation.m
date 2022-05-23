clear
clc
clear global
close all

[rgb,rgb1] = mysettings;

global ratetest
global c_d

%% loading synthetic experimental data
load synthetic_expt_data.mat
exptdata{1,1} = SynData{1,1};   % C/2
exptdata{2,1} = SynData{2,1};   % 1C
exptdata{3,1} = SynData{3,1};   % 2C
exptdata{4,1} = SynData{4,1};   % 3C

%% Comparing p2D model simulation data (performed using SUNDIALS-IDA)
load initialguess_p2D_modelsim.mat
p2d_t = InitialGuess_1C_Dch(:,1);
p2d_v = InitialGuess_1C_Dch(:,2);


%%  List of parameters in Tank-in-Series Battery Model
p = LiTank_Parameters;
capacity = p(24);
p(1) = exp(p(1));
p(2) = exp(p(2));
ratetest = -1;
Actual_Capacity = 1.78; %Ah

%% Initial Conditions
y0 = LiTank_ICs(p);

%% Solving using ode15s
% MASS MATRIX
odevars = 6;           % Number of ODE variables
aevars = 8;            % Number of AE variables
node = odevars; nae = aevars;
M = [eye(node),zeros(node,nae);zeros(nae,node+nae)];
M = sparse(M);
Absolute_Tol = 1e-8*ones(1,odevars+aevars);
% Options for CC (discharge/charge)
options_cc.Mass         = M;
options_cc.MassSingular = 'yes';
options_cc.Events       = @stop_condition_cc;
options_cc.RelTol       = 1E-6;
options_cc.AbsTol       = Absolute_Tol;
% options_cc.Stats        = 'on';
% options_cc.MaxStep      = 500;
% options_cc.BDF          = 'on';
% options_cc.InitialStep  = 1E-2;
% options_cc.MaxOrder     = 2;

%% Current density and final time for discharge
ft_dch = abs(3600/ratetest)*3;
tspan_dch = linspace(0,ft_dch,500);%[0,ft_dch];
iapp_dch = ratetest*capacity;

% c_d is the option to define stop conditions based on charge or discharge.
% if charging simulation is performed then set c_d = 1;
% if discharging simulation is performed then set c_d = 2;
c_d = 2;
c_rate = [-0.5,-1,-2,-3];
for i = 1:numel(c_rate)
    fprintf('Discharge Simulation : %f C-rate \n',c_rate(i))
    ratetest = c_rate(i);
    tic
    [t,y] = ode15s(@(t,y)LiTank_Dyn(t,y,p),tspan_dch,y0,options_cc);
    toc
    simdata{i,1} = [t,y,y(:,6)*(Actual_Capacity/capacity)];
end

figure('units','normalized','outerposition',[0 0 1 1])
plot(simdata{1,1}(:,16),simdata{1,1}(:,15),'Color',rgb.darkgreen)
hold on
plot(exptdata{1,1}(:,3),exptdata{1,1}(:,2),'Color',rgb.wine,'LineStyle','-.')

plot(simdata{2,1}(:,16),simdata{2,1}(:,15),'Color',rgb.darkolivegreen)
plot(exptdata{2,1}(:,3),exptdata{2,1}(:,2),'Color',rgb.crimson,'LineStyle','-.')

plot(simdata{3,1}(:,16),simdata{3,1}(:,15),'Color',rgb.navyblue)
plot(exptdata{3,1}(:,3),exptdata{3,1}(:,2),'Color',rgb.orangered,'LineStyle','-.')

plot(simdata{4,1}(:,16),simdata{4,1}(:,15),'Color',rgb.slateblue)
plot(exptdata{4,1}(:,3),exptdata{4,1}(:,2),'Color',rgb.darkgoldenrod,'LineStyle','-.')
xlabel('Time (s)')
ylabel('Cell Potential (V)')
kk = legend('TiS - C/2','Synthetic Expt Data - C/2','TiS - 1C','Synthetic Expt Data - 1C','TiS - 2C','Synthetic Expt Data - 2C','TiS - 3C','Synthetic Expt Data - 3C');
kk.EdgeColor = 'none';
pbaspect([1 1 1])

%% Plotting 1C rate TiS model internal states

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultLineMarkerSize',6)
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultTextFontSize',12);
set(0,'DefaultAxesFontWeight','bold');
set(0,'DefaultAxesLineWidth',2);

figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,3,1),plot(simdata{2,1}(:,1),simdata{2,1}(:,15),'Color',rgb.navyblue)
hold on
plot(exptdata{2,1}(:,1),exptdata{2,1}(:,2),'Color',rgb.darkgreen)
xlabel('Time (s)')
ylabel('Cell Potential (V)')
kk = legend('1 C-rate');
kk.EdgeColor = 'none';

subplot(2,3,2),plot(simdata{2,1}(:,1),simdata{2,1}(:,2),'Color',rgb.wine)
hold on
subplot(2,3,2),plot(simdata{2,1}(:,1),simdata{2,1}(:,3),'Color',rgb.darkgreen)
subplot(2,3,2),plot(simdata{2,1}(:,1),simdata{2,1}(:,4),'Color',rgb.navyblue)
xlabel('{\boldmath $Time (s)$}','Interpreter','LaTeX')
ylabel('{\boldmath $Electrolyte Conc \times 10^{3}  (mol/m^{3})$}','Interpreter','LaTeX')
kk=legend('Positive Electrode','Separator','Negative Electrode');
kk.Location = 'southwest';

subplot(2,3,3),plot(simdata{2,1}(:,1),log(simdata{2,1}(:,8))*1E3,'Color',rgb.wine)
hold on
subplot(2,3,3),plot(simdata{2,1}(:,1),log(simdata{2,1}(:,9))*1E3,'Color',rgb.darkgreen)
subplot(2,3,3),plot(simdata{2,1}(:,1),log(simdata{2,1}(:,10))*1E3,'Color',rgb.navyblue)
xlabel('{\boldmath $Time (s)$}','Interpreter','LaTeX')
ylabel('{\boldmath $Electrolyte Potential (mV)$}','Interpreter','LaTeX')

subplot(2,3,4),plot(simdata{2,1}(:,1),simdata{2,1}(:,5),'Color',rgb.wine)
hold on
subplot(2,3,4),plot(simdata{2,1}(:,1),simdata{2,1}(:,6),'Color',rgb.navyblue)
subplot(2,3,4),plot(simdata{2,1}(:,1),simdata{2,1}(:,13),'Color',rgb.wine*2)
subplot(2,3,4),plot(simdata{2,1}(:,1),simdata{2,1}(:,14),'Color',rgb.navyblue*2)
xlabel('{\boldmath $Time (s)$}','Interpreter','LaTeX')
ylabel('{\boldmath $Concentration_{Solid} (mol/m^{3})$}','Interpreter','LaTeX')
kk=legend('Positive Electrode (Avg)','Negative Electrode (Avg)',...
    'Positive Electrode (Surface)','Negative Electrode (Surface)');
kk.Location = 'southwest';
 
subplot(2,3,5),plot(simdata{2,1}(:,1),simdata{2,1}(:,11),'Color',rgb.wine)
hold on
plot(simdata{2,1}(:,1),simdata{2,1}(:,12),'Color',rgb.navyblue)
xlabel('{\boldmath $Time (s)$}','Interpreter','LaTeX')
ylabel('{\boldmath $\phi_{1} (V)$}','Interpreter','LaTeX')
legend('Positive Electrode','Negative Electrode')
 
subplot(2,3,6),plot(simdata{2,1}(:,1),simdata{2,1}(:,7)*(Actual_Capacity/capacity),'Color',rgb.navyblue)
xlabel('{\boldmath $Time (s)$}','Interpreter','LaTeX')
ylabel('{\boldmath $Q (Ah)$}','Interpreter','LaTeX')

