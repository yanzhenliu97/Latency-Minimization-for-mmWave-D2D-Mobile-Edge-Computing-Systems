%Distance Dy plot
load('performace_distance_user_BS_02.mat')
proposed_delay = proposed_delay*1e3;
heuristic_delay = heuristic_delay*1e3;
local_computing_delay = local_computing_delay*1e3;
edge_computing_delay = edge_computing_delay*1e3;
AO_delay = AO_delay*1e3;
OMP_delay = OMP_delay*1e3;
CM_delay = CM_delay*1e3;

size = 12;
figure()
box on
set(gca,'fontsize',size)
set(gca,'Xcolor','k')
xlabel('Distance between the BS and users D_y (m)','fontsize',size)
ylabel('System delay (ms)','fontsize',size)
set(gca,'XLim',[20 170]);
set(gca,'xTick',[20:30:170])
% set(gca,'YLim',[2 12]);
% set(gca,'yTick',[2:2:12]);
x = 20:30:170;
proposed_plot=line(x,mean(proposed_delay,1));
set(proposed_plot,'color','b')
set(proposed_plot,'linestyle','-')
set(proposed_plot,'linewidth',1.5)
set(proposed_plot,'marker','o')
hold on
heuristic_plot=line(x,mean(heuristic_delay,1));
set(heuristic_plot,'color','c')
set(heuristic_plot,'linestyle','-')
set(heuristic_plot,'linewidth',1.5)
set(heuristic_plot,'marker','s')
hold on
binary_plot=line(x,mean(min(local_computing_delay,edge_computing_delay),1));
set(binary_plot,'color','k')
set(binary_plot,'linestyle','-')
set(binary_plot,'linewidth',1.5)
set(binary_plot,'marker','d')
hold on
AO_plot=line(x,mean(AO_delay,1));
set(AO_plot,'color','m')
set(AO_plot,'linestyle','-')
set(AO_plot,'linewidth',1.5)
set(AO_plot,'marker','>')
hold on
CM_plot=line(x,mean(CM_delay,1));
set(CM_plot,'color',[1,0.5,0])
set(CM_plot,'linestyle','-')
set(CM_plot,'linewidth',1.5)
set(CM_plot,'marker','<')
hold on
OMP_plot=line(x,mean(OMP_delay,1));
set(OMP_plot,'color','g')
set(OMP_plot,'linestyle','-')
set(OMP_plot,'linewidth',1.5)
set(OMP_plot,'marker','^')

grid on
legend([proposed_plot,heuristic_plot,binary_plot,OMP_plot,CM_plot,AO_plot],'Two-timescale proposed','Two-timescale heuristic beamforming','Two-timescale binary offloading','Single-timescale OMP','Single-timescale CM','Single-timescale AO','fontsize',9)
