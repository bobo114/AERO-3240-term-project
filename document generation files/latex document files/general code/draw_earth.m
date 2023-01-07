%The following code must be run in matlab live after calculating r_I in
%PROBA2.m

figure

grid on

earth_sphere(gca,'km')

set(gca,'FontSize',9,'FontName', 'Times')
title('PROBA 2 Orbital Path in ECI Reference Frame')
xlabel('rx_I (km)')
ylabel('ry_I (km)')
zlabel('rz_I (km)')
axis equal
view([62 32])
for i = 1:500:length(r_I(:,1))
    hold on
    plot3(r_I(1:i,1),r_I(1:i,2),r_I(1:i,3),'k','linewidth', 2)
    pause(0.00001);
    drawnow;
end
plot3(r_I(:,1),r_I(:,2),r_I(:,3),'k','linewidth', 2)
drawnow;