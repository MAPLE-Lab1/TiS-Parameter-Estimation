clear
clc
clear global
close all
set(0,'DefaultLineLineWidth',3)
set(0,'DefaultLineMarkerSize',6)
set(0,'DefaultAxesFontSize',20);
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultTextFontSize',20);
set(0,'DefaultAxesFontWeight','bold');
set(0,'DefaultAxesLineWidth',2);
warning off

tic
global exptdata 

%% Estimating Transport Parameters in scaled form
% Loading synthetic experimental data
load synthetic_expt_data.mat
exptdata{1,1} = SynData{1,1}; % C/2
exptdata{2,1} = SynData{2,1}; % 1C
exptdata{3,1} = SynData{3,1}; % 2C
exptdata{4,1} = SynData{4,1}; % 3C

%% Defining Bounds
%             LB            UB
Bounds = [log(1E-15)  log(9E-13);   % Dsn
          log(1E-15)  log(9E-13)];  % Dsp     
LB = Bounds(:,1)'; % Lower Bounds
UB = Bounds(:,2)'; % Upper Bounds
pars_names = {'Dsn_optimal','Dsp_optimal'};

%% Initial guess to optimizer
p0 = [log(1.4E-13);  % Dsn
      log(2E-13)];   % Dsp
  
%% Options to optimizer
options = optimoptions('ga','MutationFcn',@mutationadaptfeasible,...
    'display','iter',...
    'SelectionFcn',@selectionremainder,...
    'Vectorized','off',...
    'Generations',50,...
    'PopulationSize',100);
% There are no linear constraints in this demo.
A1 = [];
B1 = []; 

%% Case - 1 : Estimationg Dsp, Dsn using C/2 & 1C data
% Model Cell potential using initial guess
Npars = numel(p0);
obj_guess = round(TankModel_ObjFcn(p0))
legname1 = strcat('$Guess, RMSE =',mat2str(obj_guess),'mV$');
load('modelsim.mat')
modelsim_guess = modelsim;
rmse_guess = rmse;
% Performing optimization to estimate parameters
[p_opt,FVAL,EXITFLAG,output,population,scores] = ga(@(LiPars) TankModel_ObjFcn(LiPars),Npars,...
    A1,B1,[],[],LB,UB,[],options);
p_estimated = exp(p_opt);


obj_conv = round(TankModel_ObjFcn(p_opt));
legname2 = strcat('$Optimal, RMSE =',mat2str(obj_conv),'mV$');
load modelsim.mat
modelsim_opt = modelsim;
rmse_opt = rmse;

for i=1:2
    fprintf('%s : %e m^/s \n',pars_names{i},p_estimated(i))
end

titlename = 'C/2';
Estimation_Results_Plots(exptdata{1,1},rmse_guess(1),modelsim_guess{1,1},rmse_opt(1),modelsim_opt{1,1},titlename)

titlename = '1C';
Estimation_Results_Plots(exptdata{2,1},rmse_guess(2),modelsim_guess{1,2},rmse_opt(2),modelsim_opt{1,2},titlename)

titlename = '2C';
Estimation_Results_Plots(exptdata{3,1},rmse_guess(3),modelsim_guess{1,3},rmse_opt(3),modelsim_opt{1,3},titlename)

titlename = '3C';
Estimation_Results_Plots(exptdata{4,1},rmse_guess(4),modelsim_guess{1,4},rmse_opt(4),modelsim_opt{1,4},titlename)



