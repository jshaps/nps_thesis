%animation for rotation plot showing reeling in tether

%first run: shapiro_thesis_main.m
% next run 
 close all
 clc
 clearvars test_drv_rot

odeopts=odeset('RelTol',1E-8,'AbsTol',1E-9);
omega = reshape(w_sim,3,length(w_sim))';
    C0=eye(2,2);
    told=tout;
    [t, C]=ode45(@rotationMatrix,tout,C0,odeopts,tout,omega);
    
% build rotation matrix at each time step from integration of w 
i = 0;
    for k = 1:100:length(C)
        i=i+1;
        theta(k) = C(k,2);
        R(:,:,i) = [cos(theta(k)) -sin(theta(k)); sin(theta(k)) cos(theta(k))];
    end
    
    % oscillation manuever
   t_rot=VideoWriter('tether_rotate.avi');
   t_rot.FrameRate = 10;
   open(t_rot)
   
   % oscillation graph manuever
   t_osc=VideoWriter('oscillate.avi');
   t_osc.FrameRate = 10;
   open(t_osc)
    
%figure
k=0;
for i = 1:100:length(C)
    k=k+1;
    
        figure(1)
        clf
        drv_centerX = 0;
        drv_centerY = L_drv(i);
        test_drv_rot(:,k) = R(:,:,k)*[0;L_drv(i)];
        drv_pos = [drv_centerX;drv_centerY];
        drv_rot = R(:,:,k)*drv_pos;
            scatter(drv_rot(1,:),drv_rot(2,:),5,'b','filled') % plot drv position

            hold on    
            plot([0,drv_rot(1,:)],[0,drv_rot(2,:)]);        % plot tether length
            scatter(0,0,5,'k','filled')     % plot center of mass position   

        hst_centerX = 0;
        hst_centerY = -L_hst(i);
        hst_pos = [hst_centerX;hst_centerY];
        hst_rot = R(:,:,k)*hst_pos;
            scatter(hst_rot(1,:),hst_rot(2,:),10,'r','filled') % plot drv position
            plot([0,hst_rot(1,:)],[0,hst_rot(2,:)]);        % plot tether length
       axis([-50 50 -10 50])
       xlabel('Departure from LV orientation')      %what are these units?
       xlabel('Departure from LV orientation')      %what are these units?
       drawnow
            f = getframe(gcf);
            writeVideo(t_rot,f);
   %xlabel('Tether Length [km]')
   %ylabel('Tether Length [km]')
   

   figure(2)
        plot(test_drv_rot(1,:));
        xlabel('time steps')
        axis([0 300 -50 50]) 
        drawnow
            g = getframe(gcf);
            writeVideo(t_osc,g);

end
    close(t_rot)
    close(t_osc)


%-----------process magnitudes of things
H = reshape(H_sim,3,length(H_sim));

    for a1=1:numel(tout)
        omegamag(a1)=norm(omega(a1,:));
        H_calc(a1,:)=(J_sys(:,:,a1)*omega(a1,:)')';
        Hmag(a1) = norm(H(:,a1));
        H_calc_mag(a1) = norm(H_calc(a1,:));
  %      KE(a1)=1/2*dot(omega(a1,:),J_sys(:,:,a1)*omega(a1,:)');
    end

% show omega, H, tether length 
% non animated plot processing
figure();
h(1) = subplot(2,2,1);
    plot(tout,L_drv,tout,L_hst,'Linewidth',2);
    ylabel('Tether Length to COM [km]')
    xlabel('Time [s]')
    grid on
h(2) = subplot(2,2,2);
    plot(tout,omegamag,'Linewidth',2,'color','g')
    ylabel('\omega magnitude [rad/s]')
    xlabel('Time [s]')
    grid on
h(3) = subplot(2,2,3); % the last (odd) axes
pos = get(h,'Position');
new = mean(cellfun(@(v)v(1),pos(1:2)));
set(h(3),'Position',[new,pos{end}(2:end)])
    plot(tout,Hmag,'Linewidth',2)
    hold on
    plot(tout,H_calc_mag,'Linewidth',2)
    ylabel('Momentum magnitude [kg rad/s^2]')
    xlabel('Time [s]')
    legend('H output from simulation','calculated H = J\omega','Location','best')
    grid on
   
function dC=rotationMatrix(t,C,timevec,Omega)
    
    omega=interp1(timevec,Omega,t);
    %Actually only need three parameters
%    dC=zeros(9,1);
    C=reshape(C,[2 2]);
    A=skewMatrix(omega)*C;
    dC=A(:);
    
end

function A=skewMatrix(omega)

    A=[0 -omega(3);omega(3) 0];
    
end    