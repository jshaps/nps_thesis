% develop animation and graphs for rotation solution


close all
clc

%SetGraphics(); 

% tether reel in manuever
   reel_vid=VideoWriter('reel_vid.avi');
   reel_vid.FrameRate = 10;
   open(reel_vid)
figure()

for i = 1: 250:length(L_drv)
    clf
    drv_centerX = 0;
    drv_centerY = L_drv(i);
    scatter(drv_centerX,drv_centerY,5,'b','filled') % plot drv position
    hold on
    scatter(0,0,5,'k','filled')     % plot center of mass position
    xlim([-.1 .1])
    ylim([-5 50])
    ylabel('Tether Length [km]')
    line([0 0],[0 L_drv(i)],'Color','k');
    line([0 0],[0 -L_hst(i)],'Color','k','LineStyle','--');
   
    hst_centerX = 0;
    hst_centerY = -L_hst(i);
    scatter(hst_centerX, hst_centerY, 10,'r','filled'); % plot hst position
%    legend('debris vehicle','COM','debris vehicle tether','debris tether','debris')
    drawnow 
        f = getframe(gcf);
        writeVideo(reel_vid,f);
    
end
    close(reel_vid)

    

