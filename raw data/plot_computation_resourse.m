%edge computation resourse plot
load('performance_computation_resourse_100mW.mat')
proposed_delay_100= proposed_delay*1e3;
heuristic_delay_100 = heuristic_delay*1e3;
binary_delay_100 = binary_delay*1e3;
AO_delay_100 = AO_delay*1e3;
OMP_delay_100 = OMP_delay*1e3;
load('performance_computation_resourse_200mW.mat')
proposed_delay_200= proposed_delay*1e3;
heuristic_delay_200 = heuristic_delay*1e3;
binary_delay_200 = binary_delay*1e3;
AO_delay_200 = AO_delay*1e3;
OMP_delay_200 = OMP_delay*1e3;
size = 12;
figure()
box on
set(gca,'fontsize',size)
set(gca,'Xcolor','k')
xlabel('computation resourse ratio \eta (FL = 200MHz)','fontsize',size)
ylabel('System delay (ms)','fontsize',size)

x = 2:2:12;
proposed100_plot=line(x,mean(proposed_delay_100,1));
set(proposed100_plot,'color','b')
set(proposed100_plot,'linestyle','-')
set(proposed100_plot,'linewidth',1.5)
set(proposed100_plot,'marker','o')
% hold on
% heuristic_plot=line(x,mean(OMP_delay_100,1));
% set(heuristic_plot,'color','g')
% set(heuristic_plot,'linestyle','-')
% set(heuristic_plot,'linewidth',1.5)
% set(heuristic_plot,'marker','s')
hold on
binary100_plot=line(x,mean(binary_delay_100,1));
set(binary100_plot,'color','k')
set(binary100_plot,'linestyle','-')
set(binary100_plot,'linewidth',1.5)
set(binary100_plot,'marker','d')
hold on
proposed200_plot=line(x,mean(proposed_delay_200,1));
set(proposed200_plot,'color','b')
set(proposed200_plot,'linestyle','-.')
set(proposed200_plot,'linewidth',1.5)
set(proposed200_plot,'marker','o')
% hold on
% CM_plot=line(x,mean(OMP_delay_200,1));
% set(CM_plot,'color','g')
% set(CM_plot,'linestyle','-.')
% set(CM_plot,'linewidth',1.5)
% set(CM_plot,'marker','s')
hold on
binary200_plot=line(x,mean(binary_delay_200,1));
set(binary200_plot,'color','k')
set(binary200_plot,'linestyle','-.')
set(binary200_plot,'linewidth',1.5)
set(binary200_plot,'marker','d')

grid on
legend([proposed100_plot,proposed200_plot,binary100_plot,binary200_plot],'Two-timescale proposed PUA=100mW','Two-timescale proposed PUA=200mW','Two-timescale binary offloading PUA=100mW','Two-timescale binary offloading PUA=200mW','fontsize',9)
