function obj = calc_rmse(texp,vexp,tmodel,vmodel)

err1 = [];
tfinal_expt = texp(end);
tfinal_model= tmodel(end);
vmodel_inter = interp1(tmodel,vmodel,texp,'linear','extrap');
err = (vmodel_inter - vexp);
% if (tfinal_model >= tfinal_expt)
%     %         display('Simulation time is greater than experimential time')
%     vmodel_inter = interp1(tmodel,vmodel,texp,'linear','extrap');
%     err = (vmodel_inter - vexp);
% elseif (tfinal_model < tfinal_expt)
%     %         display('Experimental time is greater than simulation time')
%     vexp_inter = interp1(texp,vexp,tmodel,'linear','extrap');
%     err = (vmodel - vexp_inter);
% end
err1 = [err1;err];
nn = length(err1);
obj = sqrt(sum(real(err1(2:end).^2))/nn);


end