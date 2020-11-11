%用于统计用户间距离d3对通信容量的影响
%矩阵维度参数
Nu=32;
Na=8
Nb=8;
Nrf1=4;
Nrf2=4;

%L是路径数
L=15; 
%噪声方差
sigma=10^(-6); 
%功率p
p=1;
Pmax=1000;
Ni=100;

H1_cell=cell(Ni,6);
H2_cell=cell(Ni,6);
H3_cell=cell(Ni,6);
downlinkrate=zeros(Ni,6);
uplinkrate=zeros(Ni,6);
d2dlinkrate=zeros(Ni,6);
CMdownlinkrate=zeros(Ni,6);
CMuplinkrate=zeros(Ni,6);
CMd2dlinkrate=zeros(Ni,6);
%计算下行信道的最大速率


for d=10:10:60
    counter=d/10
%     if Nu==40
%         start_point=30;
%     else
%         start_point=1;
%     end
i=1;
    while 1 
        %产生信道
        
H1=mmWavechannel(Nu,Na,L,120);
H2=mmWavechannel(Nb,Nu,L,120);
H3=mmWavechannel(Nb,Na,L,d);
H1_cell{i,counter}=H1;
H2_cell{i,counter}=H2;
H3_cell{i,counter}=H3;

[uplinkrate(i,counter),CMuplinkrate(i,counter)]=pcccp_uplink(Nu,Na,Nrf1,H1,sigma);
[downlinkrate(i,counter),CMdownlinkrate(i,counter)]=pcccp_downlink(Nu,Nb,Nrf2,H2,Pmax,sigma);
[d2dlinkrate(i,counter),CMd2dlinkrate(i,counter)]=D2Dlink(Na,Nb,H3,sigma);

if(CMd2dlinkrate(i,counter)<4)
    i=i-1;  
end
i=i+1;
if(i>=Ni)
    break;
end
    end
end