function [value,isterminal,direction] = stop_condition_cc(t,y)

global c_d

if(c_d==1)
    value = y(14)-4.2;  % charging at CC 
elseif(c_d==2)
    value = y(14)-2.8;  % discharging at CC 
end
isterminal = 1; % stop the integration
direction = 0;  % negative direction
return;