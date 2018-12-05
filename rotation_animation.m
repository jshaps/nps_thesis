%animation for rotation plot showing reeling in tether

%first run: shapiro_thesis_main.m
% next run 
close all
clc

odeopts=odeset('RelTol',1E-8,'AbsTol',1E-9);
omega = reshape(w_sim,3,length(w_sim))';
    C0=eye(2,2);
    told=tout;
    [t, C]=ode45(@rotationMatrix,tout,C0,odeopts,tout,omega);
    
% build rotation matrix at each time step from integration of w    
    for k = 1:length(C)
        theta(k) = C(k,2);
        R(:,:,k) = [cos(theta(k)) -sin(theta(k)); sin(theta(k)) cos(theta(k))];
    end
    
figure
for i = 1:200:length(C)
    clf
   x = cos(C(i,3));
   y = sin(C(i,3));
    drv_centerX = L_drv(i);
    drv_centerY = L_drv(i);
    drv_pos = [drv_centerX;drv_centerY];
            drv_rot = R(:,:,i)*drv_pos;
            test_drv_rot(:,i) = drv_rot;
        scatter(drv_pos(1,:),drv_pos(2,:),5,'b','filled')
        scatter(drv_rot(1,:),drv_rot(2,:),5,'b','filled') % plot drv position
    hold on    
        plot([0,drv_rot(1,:)],[0,drv_rot(2,:)]);        % plot tether length
    
        scatter(0,0,5,'k','filled')     % plot center of mass position   
        
    hst_centerX = -L_hst(i);
    hst_centerY = -L_hst(i);
        hst_pos = [hst_centerX;hst_centerY];
        
        hst_rot = R(:,:,i)*hst_pos;
                scatter(hst_rot(1,:),hst_rot(2,:),10,'r','filled') % plot drv position
        plot([0,hst_rot(1,:)],[0,hst_rot(2,:)]);        % plot tether length
    

   axis([-20 70 -20 70])
   xlabel('Tether Length [km]')
   ylabel('Tether Length [km]')
   drawnow

end
    
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