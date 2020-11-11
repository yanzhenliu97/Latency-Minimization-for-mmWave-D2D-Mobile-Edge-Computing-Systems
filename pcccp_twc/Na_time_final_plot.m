load Na_result.mat

W=100e6; %带宽
L=10e6; %要计算的比特数
Fe=6000e6; %边缘计算服务器速率
FL=1800e6; %本地计算服务器速率
alpha=0.1; %计算压缩率

x=2:2:12;
for i=1:6
    for j=1:50 
k1=L/(W*uplinkrate(j,i));
ke=L/Fe;
k2=alpha*L/(W*downlinkrate(j,i));
kL=L/FL;
k3=alpha*L/(W*d2dlinkrate(j,i));
[rho(j,i),optime(j,i),condition(j,i)]=best_partition(k1,ke,k2,kL,k3);
localtime(j,i)=kL+k3;
edgetime(j,i)=k1+ke+k2;

%计算cm算法的时间
k1=L/(W*CMuplinkrate(j,i));
k2=alpha*L/(W*CMdownlinkrate(j,i));
k3=alpha*L/(W*CMd2dlinkrate(j,1));
%[cmrho(j,i),cmtime(j,i),cmcondition(j,i)]=best_partition(k1,ke,k2,kL,k3);
[cmtime2(j,i)]=rho_partition(k1,ke,k2,kL,k3,0.5);
    end
end
optimal_time=1000*mean(optime);
optimal_rho=mean(rho);
localcomputing_time=1000*mean(localtime);
edgecomputing_time=1000*mean(edgetime);
%cm_optimal_time=1000*mean(cmtime);
cm_optimal_time2=1000*mean(cmtime2);
size=14;
figure()
%[ax,time_plot,rho_plot]=plotyy(x,optimal_time(1:6),x,optimal_rho(1:6));

% set(gca,'fontsize',size-3)
% set(gca,'Xcolor','k')
% set(gca,'XLim',[250 1500]);
% set(gca,'xTick',[250:250:1500])

% set(ax(1),'YLim',[5 35]);
% set(ax(1),'yTick',[5:5:35]);
% set(ax(1),'Ycolor','k');
% set(ax(2),'fontsize',size-3);
% set(ax(2),'Ycolor','k');
% set(get(ax(1),'Ylabel'),'String','System delay (ms)','fontsize',size) 
% set(get(ax(2),'Ylabel'),'String','Optimal offloading ratio','fontsize',size)
% xlabel('F_E (Mbps)','fontsize',size-1)
set(gca,'fontsize',size-3)
set(gca,'Xcolor','k')
xlabel('N_a','fontsize',size-1)
ylabel('System delay (ms)','fontsize',size)
set(gca,'XLim',[2 12]);
set(gca,'xTick',[2:2:12])
set(gca,'YLim',[3 10]);
set(gca,'yTick',[3:2:11]);

time_plot=line(x,optimal_time(1:6));
set(time_plot,'color','b')
set(time_plot,'linestyle','-')
set(time_plot,'linewidth',2)
set(time_plot,'marker','o')
% set(rho_plot,'color','r')
% set(rho_plot,'linestyle','none')
% set(rho_plot,'linewidth',2)
% set(rho_plot,'marker','p')
hold on
local_plot=line(x,localcomputing_time(1:6));
set(local_plot,'color','c')
set(local_plot,'linestyle','-')
set(local_plot,'linewidth',2)
set(local_plot,'marker','s')
hold on
edge_plot=line(x,edgecomputing_time(1:6));
set(edge_plot,'color','k')
set(edge_plot,'linestyle','-')
set(edge_plot,'linewidth',2)
set(edge_plot,'marker','d')
hold on
cm2_plot=line(x,cm_optimal_time2(1:6));
set(cm2_plot,'color','m')
set(cm2_plot,'linestyle','-')
set(cm2_plot,'linewidth',2)
set(cm2_plot,'marker','>')
% cm_plot=line(x,cm_optimal_time(1:6));
% set(cm_plot,'color','r')
% set(cm_plot,'linestyle','-')
% set(cm_plot,'linewidth',2)
% set(cm_plot,'marker','>')
grid on
legend([time_plot,local_plot,edge_plot,cm2_plot],'Proposed joint design algorithm','Local computing','Edge computing','Channel matching \rho=0.5','fontsize',size-2.5)
