function [rgbcolor,rgbcolor1] = mysettings
set(0,'DefaultLineLineWidth',2)                 % Change this number if the plot line looks thick or thin
set(0,'DefaultLineMarkerSize',6)               % Size of the symbols that are used in the plots
set(0,'DefaultAxesFontSize',20);
% set(0,'DefaultAxesFontName', 'Times New Roman') % Use either Times New Roman or Arial
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultTextFontSize',24);                % Font of the axes labels
set(0,'DefaultAxesFontWeight','bold');          % Option to bold axes labels
set(0,'DefaultAxesLineWidth',4);

rgbcolor.wine = [0.474509803921569 0 0.0823529411764706];    % 1
rgbcolor.maroon = [128,0,0]./255;                            % 2
rgbcolor.orangered = [255,69,0]./255;                        % 3
rgbcolor.olive = [128,128,0]./255;                           % 4
rgbcolor.darkolivegreen = [85,107,47]./255;                  % 5
rgbcolor.seagreen = [46,139,87]./255;                        % 6
rgbcolor.midnightblue = [25,25,112]./255;                    % 7
rgbcolor.indigo = [75,0,130]./255;                           % 8
rgbcolor.purple = [128,0,128]./255;                          % 9
rgbcolor.saddlebrown = [139,69,19]./255;                     % 10
rgbcolor.darkgreen = [0 0.349019607843137 0];                % 11
rgbcolor.navyblue = [0 0.0980392156862745 0.482352941176471];% 12
rgbcolor.darkcyan = [0.00,0.55,0.55];                        % 13
rgbcolor.peach = [251 111 66] ./ 255;                        % 14
rgbcolor.teal = [18 150 155] ./ 255;                         % 15
rgbcolor.darkgoldenrod = [184,134,11]./255;                  % 16
rgbcolor.dark_turquoise = [0.00,0.81,0.82];                  % 17
rgbcolor.dark_violet =[0.58,0.00,0.83];                      % 18
rgbcolor.coral = [255,127,80]./255;                          % 19
rgbcolor.darkorange = [255,140,0]./255;                      % 20
rgbcolor.crimson = [220,20,60]./255;                         % 21
rgbcolor.darkslategray = [47,79,79]./255;                    % 22
rgbcolor.slateblue = [106,90,205]./255;                      % 23
rgbcolor.royalblue = [65,105,225]./255;                      % 24
%**********************************************************%


rgbcolor1(1,:) = rgbcolor.darkgreen;
rgbcolor1(2,:) = rgbcolor.navyblue;
rgbcolor1(3,:) = rgbcolor.wine;

rgbcolor1(4,:) = rgbcolor.darkolivegreen;
rgbcolor1(5,:) = rgbcolor.midnightblue;
rgbcolor1(6,:) = rgbcolor.orangered;

rgbcolor1(7,:) = rgbcolor.olive;
rgbcolor1(8,:) = rgbcolor.indigo;
rgbcolor1(9,:) = rgbcolor.peach;

rgbcolor1(10,:) = rgbcolor.saddlebrown;
rgbcolor1(11,:) = rgbcolor.dark_turquoise;
rgbcolor1(12,:) = rgbcolor.coral;

rgbcolor1(13,:) = rgbcolor.darkslategray;
rgbcolor1(14,:) = rgbcolor.slateblue;
rgbcolor1(15,:) = rgbcolor.teal;

rgbcolor1(16,:) = 2*rgbcolor.darkgreen;
rgbcolor1(17,:) = 2*rgbcolor.navyblue;
rgbcolor1(18,:) = 2*rgbcolor.wine;

rgbcolor1(19,:) = 2*rgbcolor.darkolivegreen;
rgbcolor1(20,:) = 2*rgbcolor.midnightblue;
rgbcolor1(21,:) = 2*rgbcolor.orangered;

end 