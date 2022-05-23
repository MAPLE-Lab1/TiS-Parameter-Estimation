function obj = TankModel_ObjFcn(LiPars)

global exptdata c_d ratetest

error = [];

%%  List of parameters in Tank-in-Series Battery Model
p = LiTank_Parameters;
p(1) = exp(LiPars(1));
p(2) = exp(LiPars(2));

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

%% Current density and final time for discharge
% c_d is the option to define stop conditions based on charge or discharge.
% if charging simulation is performed then set c_d = 1;
% if discharging simulation is performed then set c_d = 2;
c_rate = [-0.5,-1,-2,-3];
try
    c_d = 2; 
    for i = 1:numel(c_rate)
        ratetest = c_rate(i);
        ft_dch = abs(3600/ratetest)+1000;
        tspan_dch = [0 ft_dch];
        [t,y] = ode15s(@(t,y)LiTank_Dyn(t,y,p),tspan_dch,y0,options_cc);
        tmodel = t;
        vmodel = y(:,14);
        rmse(i) = calc_rmse(exptdata{i,1}(:,1),exptdata{i,1}(:,2),tmodel,vmodel)*1000;
        modelsim{i}.t = t;
        modelsim{i}.pot = y(:,14);
        modelsim{i}.statevars = y;
        save modelsim modelsim rmse
    end
    obj = sum(rmse);
catch
        obj = 1E4;
        modelsim{1}.t = 0;
        modelsim{1}.pot = 0;

        modelsim{2}.t = 0;
        modelsim{2}.pot = 0;
        
        rmse = 0;
        save modelsim modelsim rmse

end


end