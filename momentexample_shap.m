function momentexample_shap(J_sys,t_sim)

% code by Dr. Richard B. Choroszucha
% University of Michigan, 2011
% Adapted by Jessica Shapiro

    J=J_sys;
    Jinv=J^(-1);
    M=zeros(3,1);       % must be 0 for angular momentum to be conserved
    odeopts=odeset('RelTol',1E-8,'AbsTol',1E-9);
    
    %t0=[0 200];
    t0=linspace(0,t_sim,500);%,round(t_sim)*10);
    %omega0=[0.1 -0.1 0.4];
    
    % playing with rotation vector here
    % w_bn angular rotation
    omega0 = [-3.59215606171920e-05,0.00182535310158493,0.00110242490610287];
    %omega0 = -.1*[0 0 1];
    
    [t omega]=ode45(@eulerequations,t0,omega0,odeopts,J,Jinv,M);

    H=zeros(numel(t),3);
    for a1=1:numel(t)
        omegamag(a1)=norm(omega(a1,:),2);
        H(a1,:)=(J*omega(a1,:)')';
        Hmag(a1)=norm(H(a1,:),2);
        KE(a1)=1/2*dot(omega(a1,:),J*omega(a1,:)');
    end

    figure(1)
        subplot(2,2,1)
            plot(t,omega(:,1),'k');
            xlabel('Time [s]')
            ylabel('Angular Velocity X [rad/s]')
        subplot(2,2,2)
            plot(t,omega(:,2),'k');
            xlabel('Time [s]')
            ylabel('Angular Velocity Y [rad/s]')
        subplot(2,2,3)
            plot(t,omega(:,3),'k');
            xlabel('Time [s]')
            ylabel('Angular Velocity Z [rad/s]')
        subplot(2,2,4)
            plot(t,omegamag,'k');
            xlabel('Time [s]')
            ylabel('Angular Velocity Magnitude [rad/s]')
    figure(2)
        subplot(2,2,1)
            plot(t,H(:,1),'k');
            xlabel('Time [s]')
            ylabel('Angular Momentum X [rad/s]')
        subplot(2,2,2)
            plot(t,H(:,2),'k');
            xlabel('Time [s]')
            ylabel('Angular Momentum Y [rad/s]')
        subplot(2,2,3)
            plot(t,H(:,3),'k');
            xlabel('Time [s]')
            ylabel('Angular Momentum Z [rad/s]')
        subplot(2,2,4)
            plot(t,Hmag,'k');
            xlabel('Time [s]')
            ylabel('Angular Momentum Magnitude [rad/s]')
    
    figure(3)
        plot(t,KE,'k')
        xlabel('Time [s]')
        ylabel('Kinetic Energy [J]')
    %{
    figure(4)
        Z=zeros(numel(t),3);
        
        quiver3(Z(:,1),Z(:,2),Z(:,3),omega(:,1),omega(:,2),omega(:,3),'b');
        %{
        for a1=1:numel(t)
            quiver3(Z(1:a1,1),Z(1:a1,2),Z(1:a1,3),omega(1:a1,1),omega(1:a1,2),omega(1:a1,3),'b');
            axis(0.5*[-1 1 -1 1 -1 1]);
            title(sprintf('Time = %f',t(a1)));
            view([45 90])
            drawnow;
        end
        %}
    %}
        
    C0=eye(3,3);
    told=t;
    [t C]=ode45(@rotationMatrix,t0,C0(:),odeopts,t,omega);
    
    omega=interp1(told,omega,t);
    Z=zeros(numel(t),3);
    
    numsphere=50;
    
    [S1X S1Y S1Z]=sphere(numsphere);
    [S2X S2Y S2Z]=sphere(numsphere);
    [C1X C1Y C1Z]=cylinder(0.5,numsphere);
    
    
    S1X=2*S1X;
    S1Y=2*S1Y;
    S1Z=2*S1Z;
    S2X=2*S2X;
    S2Y=2*S2Y;
    S2Z=2*S2Z;
    
    S1Y=S1Y+5;
    S2Y=S2Y-5;
    
% ,'Uncompressed AVI'    
   %AVI=VideoWriter('init_tether_50km.avi');
   AVI=VideoWriter('final_tether_10km.avi');
   AVI.FrameRate = 10;
   open(AVI)
    figure(5)
        h1=surf(S1X,S1Y,S1Z);
        hold on;
        h2=surf(S2X,S2Y,S2Z);
        inter=C1Y;
        C1Y=10*C1Z-5;
        C1Z=inter;
        h3=surf(C1X,C1Y,C1Z);
        S1Colour=1/2*get(h1,'CData');
        S2Colour=get(h2,'CData');
        C1Colour=get(h3,'CData');
    close(5);
        
    numsphere=numsphere+1;
    %,'fps',20

    
    for o1=1%:2
    
%         if (o1==2)
%             close(6);
%             h=figure(6);
%                 set(h,'Color',[1 1 1]);
%             text(0,0,0,{'What if we were to view it from','the perspective of the angular','velocity vector?'}) 
%             axis([-0.5 1 -1 1])
%             axis off;
%             for a1=1:200
%                 %F=getframe(h);
%                 %AVI=addframe(AVI,F);
%             end
%             axis on;
%         end
        
        m=10;
        centre=zeros(floor(numel(t)/m),3);
        centre0=[0 5 0]';
        
        for a1=1:m:numel(t)

            R=reshape(C(a1,:),[3 3]);
            frame=2*R'*C0;
            h=figure(6);
                set(h,'Color',[1 1 1]);

                s1=R*[S1X(:)';S1Y(:)';S1Z(:)'];
                s2=R*[S2X(:)';S2Y(:)';S2Z(:)'];
                cy1=R*[C1X(:)';C1Y(:)';C1Z(:)'];
                S1RX=reshape(s1(1,:),[numsphere numsphere]);
                S1RY=reshape(s1(2,:),[numsphere numsphere]);
                S1RZ=reshape(s1(3,:),[numsphere numsphere]);
                S2RX=reshape(s2(1,:),[numsphere numsphere]);
                S2RY=reshape(s2(2,:),[numsphere numsphere]);
                S2RZ=reshape(s2(3,:),[numsphere numsphere]);
                C1RX=reshape(cy1(1,:),[2 numsphere]);
                C1RY=reshape(cy1(2,:),[2 numsphere]);
                C1RZ=reshape(cy1(3,:),[2 numsphere]);

                centre((a1-1)/10+1,:)=(R*centre0)';
                
                surf(S1RX,S1RY,S1RZ,'LineStyle','None','CData',S1Colour);
                hold on;
                surf(S2RX,S2RY,S2RZ,'LineStyle','None','CData',S2Colour);
                surf(C1RX,C1RY,C1RZ,'LineStyle','None','CData',C1Colour);
                for a2=1:3
                    quiver3(0,0,0,frame(a2,1),frame(a2,2),frame(a2,3),'k')
                end
                plot3(10*omega(1:a1,1),10*omega(1:a1,2),10*omega(1:a1,3),'b','LineWidth',2)
                if (o1==1)
                    plot3(centre(1:((a1-1)/10+1),1),centre(1:((a1-1)/10+1),2),centre(1:((a1-1)/10+1),3),'k','LineWidth',0.5)
                end
                quiver3(Z(a1,1),Z(a1,2),Z(a1,3),omega(a1,1),omega(a1,2),omega(a1,3),10,'b');

                v=[omega(a1,1),omega(a1,2),omega(a1,3)];
                v=v/norm(v,2);

                title(sprintf('Time = %f',t(a1)));
                axis square;
                
                axis(8*[-1 1 -1 1 -1 1]);
                
                if (o1==2)
                    view(v);
                end
                %f = getframe(gcf);
                %writeVideo(AVI,f);
                drawnow;
                f = getframe(gcf);
                writeVideo(AVI,f);
                hold off
                
        end
        

    end
    
    framecount = AVI.FrameCount
    close(AVI)
    fprintf('Fin.\n');
        
end


function dx=eulerequations(t,x,J,Jinv,M)

    dx=zeros(3,1);
    dx(1:3)= -Jinv*(cross(x(1:3),J*x(1:3))+ M);
    
end


function dC=rotationMatrix(t,C,timevec,Omega)
    
    omega=interp1(timevec,Omega,t);
    %Actually only need three parameters
    dC=zeros(9,1);
    C=reshape(C,[3 3]);
    A=crossMatrix(omega)*C;
    dC=A(:);
end

function A=crossMatrix(omega)

    A=[0 -omega(3) omega(2);omega(3) 0 -omega(1);-omega(2) omega(1) 0];
end
