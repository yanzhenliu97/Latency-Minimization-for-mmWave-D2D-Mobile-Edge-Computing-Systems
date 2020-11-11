load d1d2result2.mat

W=100e6; %带宽
L=10e6; %要计算的比特数
Fe=6400e6; %边缘计算服务器速率
FL=900e6; %本地计算服务器速率
alpha=0.1; %计算压缩率

x=1600:1600:1600+1600*5;
for i=1:6
    
    Fe=x(i)*1e6;
    for j=1:50

k1=L/(W*uplinkrate(j,3));
ke=L/Fe;
k2=alpha*L/(W*downlinkrate(j,3));
kL=L/FL;
k3=alpha*L/(W*d2dlinkrate(j,3));
[rho(j,i),optime(j,i),condition(j,i)]=best_partition(k1,ke,k2,kL,k3);
localtime(j,i)=kL+k3;
edgetime(j,i)=k1+ke+k2;
[optime1(j,i)]=rho_partition(k1,ke,k2,kL,k3,0.5);
%计算cm算法的时间
k1=L/(W*CMuplinkrate(j,3));
k2=alpha*L/(W*CMdownlinkrate(j,3));
k3=alpha*L/(W*CMd2dlinkrate(j,3));
[cmrho(j,i),cmtime(j,i),cmcondition(j,i)]=best_partition(k1,ke,k2,kL,k3);
%(W*uplinkrate(j,4)+Fe)/(W*uplinkrate(j,4)+Fe+FL)
[cmtime2(j,i)]=rho_partition(k1,ke,k2,kL,k3,cmrho(j,i));
    end
end
optimal_time=1000*mean(optime);
optimal_rho=mean(rho);
localcomputing_time=1000*mean(localtime);
edgecomputing_time=1000*mean(edgetime);
%cm_optimal_time=1000*mean(cmtime);
cm_optimal_time1=1000*mean(optime1);
cm_optimal_time2=1000*mean(cmtime2);

size=14;
figure
linesize=1.5;
set(gca,'fontsize',size-3)
set(gca,'Xcolor','k')
xlabel('F_E (Mbps)','fontsize',size)
ylabel('System delay (ms)','fontsize',size)
set(gca,'XLim',[1600 9600]);
set(gca,'xTick',[1600:1600:9600])
% set(gca,'YLim',[3 13]);
% set(gca,'yTick',[3:2:13]);
size=14;
set(gca,'fontsize',size-3)
set(gca,'Xcolor','k')
xlabel('F_E (Mbps)','fontsize',size-1)
ylabel('System delay (ms)','fontsize',size)
% set(gca,'XLim',[400 1900]);
% set(gca,'xTick',[400:300:1900])
% set(gca,'YLim',[3 13]);
% set(gca,'yTick',[3:2:13]);
time_plot=line(x,optimal_time(1:6));
set(time_plot,'color','b')
set(time_plot,'linestyle','-')
set(time_plot,'linewidth',linesize)
set(time_plot,'marker','o')
% set(rho_plot,'color','r')
% set(rho_plot,'linestyle','none')
% set(rho_plot,'linewidth',2)
% set(rho_plot,'marker','p')
hold on
local_plot=line(x,localcomputing_time(1:6));
set(local_plot,'color','c')
set(local_plot,'linestyle','-')
set(local_plot,'linewidth',linesize)
set(local_plot,'marker','s')
hold on
edge_plot=line(x,edgecomputing_time(1:6));
set(edge_plot,'color','k')
set(edge_plot,'linestyle','-')
set(edge_plot,'linewidth',linesize)
set(edge_plot,'marker','d')
hold on
% cm1_plot=line(x,cm_optimal_time1(1:6));
% set(cm1_plot,'color','r')
% set(cm1_plot,'linestyle','-')
% set(cm1_plot,'linewidth',2)
% set(cm1_plot,'marker','<')
% hold on
cm2_plot=line(x,cm_optimal_time2(1:6));
set(cm2_plot,'color','m')
set(cm2_plot,'linestyle','-')
set(cm2_plot,'linewidth',linesize)
set(cm2_plot,'marker','>')
hold on
% cm_localplot=line(x,cm_localtime(1:6));
% set(cm_localplot,'color','m')
% set(cm_localplot,'linestyle','-')
% set(cm_localplot,'linewidth',2)
% set(cm_localplot,'marker','s')
% hold on
% cm_edgeplot=line(x,cm_edgetime(1:6));
% set(cm_edgeplot,'color','r')
% set(cm_edgeplot,'linestyle','-')
% set(cm_edgeplot,'linewidth',2)
% set(cm_edgeplot,'marker','s')
grid on
box on
legend([time_plot,local_plot,edge_plot,cm2_plot],'Proposed joint design algorithm','Local computing','Edge computing','Channel matching ')
set(gca,'fontsize',size-2)