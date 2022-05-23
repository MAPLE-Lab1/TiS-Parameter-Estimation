clear all
clc
close all
[rgb,~] = mysettings;

rng(0);
sigma = 0.125/100;
noise_voltage = (10/1000)*(1/10)*randi([0,10],201,1);%normrnd(0,1,[201,1]);


initial_guess = [2.82418480228790e-10,1.40000000000000e-14,2.00000000000000e-14,6.62600000000000e-10,2.40500000000000e-10,5.00000000000000e-06,5.00000000000000e-06,1200,31080,51830,0.300000000000000,0.300000000000000,0.400000000000000,0.0380000000000000,0.120000000000000,4.00000000000000e-05,2.50000000000000e-05,3.65500000000000e-05,1.50000000000000,1.50000000000000,1.50000000000000,100,100,0.359733744935366,0.790797940797941,17.5400000000000,1,2.80000000000000];

c_rate = [1/2,1,2,3];

for i = 1:numel(c_rate)
    fprintf('C-rate : %2.1f C \n',c_rate(i))
        pars = initial_guess;
        pars(27) = c_rate(i);
        
        fid = fopen('p2D_ippars.txt','w');
        fprintf(fid,'%e \n',pars);
        fclose(fid);
        [a,b]=dos('p2Dmodel.exe');

        try
            p2D_dch{i,1} = load('DischargeCurve.txt'); 
            delete('DischargeCurve.txt')
        catch
            p2D_dch{i,1} = ones(10,2);
        end
        Npts = numel(p2D_dch{i,1}(:,1));
        idx = [round(linspace(1,Npts-1,200)),Npts];
        time    = p2D_dch{i,1}(idx,1);
        voltage = p2D_dch{i,1}(idx,2)-0.005*ones(numel(idx),1)+noise_voltage;
        Q       = time*c_rate(i)*1.78/3600;        
        
        SynData{i,1} = [time,voltage,Q];
end



figure('Name', 'Cell Potential Vs Time','units','normalized','outerposition',[0 0 1 1])
plot(SynData{1,1}(:,3),SynData{1,1}(:,2),'color',rgb.wine,'LineStyle','-.')
hold on
plot(SynData{2,1}(:,3),SynData{2,1}(:,2),'color',rgb.crimson,'LineStyle','-.')
plot(SynData{3,1}(:,3),SynData{3,1}(:,2),'color',rgb.orangered,'LineStyle','-.')
plot(SynData{4,1}(:,3),SynData{4,1}(:,2),'color',rgb.darkgoldenrod,'LineStyle','-.')
xlabel('Capacity (Ah)')
ylabel('Cell Potential (V)')
kk = legend('C/2','1C','2C','3C');
kk.EdgeColor = 'none';
pbaspect([1 1 1])


save synthetic_expt_data SynData