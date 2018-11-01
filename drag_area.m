%% plots for stela data
clear 
clc
close all 

% import csv file
%area_50T = VarName7;
%life_50T = VarName15;
load('debris.mat')

sphere_r = sqrt(area*(1/pi));
sphere_r = round(sphere_r,2);
sph_label = {num2str(sphere_r(1)),num2str(sphere_r(10)),num2str(sphere_r(20)),num2str(sphere_r(30)),...
    num2str(sphere_r(40)),num2str(sphere_r(50)),num2str(sphere_r(60)),...
    num2str(sphere_r(70)),num2str(sphere_r(80)),num2str(sphere_r(90)),num2str(sphere_r(99))};


figure()
%title('Effect of increasing drag area on debris lifetime')
hAX=axes;                 % first axes, save handle
pos=get(hAX,'position');   % get the position vector
pos1=pos(2);              % save the original bottom position
pos(2)=pos(2)+pos1; pos(4)=pos(4)-pos1;  % raise bottom/reduce height->
                                         % same overall upper position
set(hAX,'position',pos)   % and resize first axes
pos(2)=pos1; pos(4)=0.01; % reset bottom to original and small height
plot(area-area(1),life,'Linewidth',2)
grid on
hold on
plot(area_50T-area(1),life_50T,'Linewidth',2)
legend('Baseline','50 km tether','Location','best')
hAX(2)=axes('position',pos,'color','none');  % and create the second
ylim(hAX(1),[0 4.5])
ylabel(hAX(1),'Lifetime [yr]')
xlabel(hAX(1),'Added Drag Area [m^2]')
xlabel(hAX(2),'Sphere radius [m]')
set(hAX(2),'xcolor','r','ycolor','r')
xticklabels(hAX(2),sph_label)


figure()
tether =(t_length(1)-t_length)*.001;
semilogy(tether,t_life)
grid on
xlabel('Tether Length [km]')
ylabel('Lifetime [yr]')
title({'Effect of tether length on debris lifetime'; '(10m sphere radius)'})