%plot_quatization
load('performance_quantization_level.mat')
proposed_delay = proposed_delay*1e3;
optimal_delay = optimal_delay*1e3;

size = 12;
figure()
box on
set(gca,'fontsize',size)
set(gca,'Xcolor','k')
xlabel('transmit power P_{UA}(mW)','fontsize',size)
ylabel('System delay (ms)','fontsize',size)
set(gca,'XLim',[50 175]);
set(gca,'xTick',[50:25:175])
% set(gca,'YLim',[2 12]);
% set(gca,'yTick',[2:2:12]);
x = 50:25:175;
proposed4_plot=line(x,mean(proposed_delay(:,:,1),1));
set(proposed4_plot,'color','g')
set(proposed4_plot,'linestyle','-')
set(proposed4_plot,'linewidth',1.5)
set(proposed4_plot,'marker','s')
hold on
optimal4_plot=line(x,mean(optimal_delay(:,:,1),1));
set(optimal4_plot,'color','g')
set(optimal4_plot,'linestyle',':')
set(optimal4_plot,'linewidth',1.5)
set(optimal4_plot,'marker','s')
hold on

proposed16_plot=line(x,mean(proposed_delay(:,:,2),1));
set(proposed16_plot,'color',[1,0.5,0])
set(proposed16_plot,'linestyle','-')
set(proposed16_plot,'linewidth',1.5)
set(proposed16_plot,'marker','*')
hold on
optimal16_plot=line(x,mean(optimal_delay(:,:,2),1));
set(optimal16_plot,'color',[1,0.5,0])
set(optimal16_plot,'linestyle',':')
set(optimal16_plot,'linewidth',1.5)
set(optimal16_plot,'marker','*')
hold on

proposed64_plot=line(x,mean(proposed_delay(:,:,3),1));
set(proposed64_plot,'color','r')
set(proposed64_plot,'linestyle','-')
set(proposed64_plot,'linewidth',1.5)
set(proposed64_plot,'marker','d')
hold on
optimal64_plot=line(x,mean(optimal_delay(:,:,3),1));
set(optimal64_plot,'color','r')
set(optimal64_plot,'linestyle',':')
set(optimal64_plot,'linewidth',1.5)
set(optimal64_plot,'marker','d')

proposedinf_plot=line(x,mean(proposed_delay(:,:,4),1));
set(proposedinf_plot,'color','b')
set(proposedinf_plot,'linestyle','-')
set(proposedinf_plot,'linewidth',1.5)
set(proposedinf_plot,'marker','o')
hold on

grid on
legend([proposed4_plot,optimal4_plot,proposed16_plot,optimal16_plot,proposed64_plot,optimal64_plot,proposedinf_plot],'Proposed (M=4)','Exhaustive search (M=4)','Proposed (M=16)','Exhaustive search (M=16)','Proposed (M=64)','Exhasuted search (M=64)','Proposed (M=inf)','fontsize',9)
