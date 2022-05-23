function Estimation_Results_Plots(exptdata,obj_guess,modelsim_guess,obj_conv,modelsim_opt,titlename)

%------Color---------------------------------------%
rgb.wine = [0.474509803921569 0 0.0823529411764706];
rgb.darkgreen = [0 0.349019607843137 0];
rgb.navyblue = [0 0.0980392156862745 0.482352941176471];
%---------------------------------------------------%

legname0 = strcat('$Expt - Synthetic Data$');
legname1 = strcat('$Guess, RMSE =',mat2str(obj_guess),'mV$');
legname2 = strcat('$Optimal, RMSE =',mat2str(obj_conv),'mV$');

figure('units','normalized','outerposition',[0 0 1 1])
plot(modelsim_guess.t/60,modelsim_guess.pot,'Color',rgb.navyblue)
hold on
title(titlename)
plot(modelsim_opt.t/60,modelsim_opt.pot,'Color',rgb.darkgreen)
plot(exptdata(:,1)/60,exptdata(:,2),'Color',rgb.wine,'LineStyle','--')
xlabel('Time (min)')
ylabel('Cell Potential (V)')
kk=legend({legname1,legname2,legname0},'Interpreter','latex');
kk.EdgeColor = 'none';
pbaspect([1 1 1])

end