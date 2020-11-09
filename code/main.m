%% coefficient initialization
FE = 1600e6  %edge computing rate
FL = 200e6  %local computing rate
L = 1e6  %bits in a time slot
alpha = 0.02  %compression ratio

B1 = 100e6  %bandwidth
B2 = 100e6
B3 = 100e6
Pmax = 30  %PA max power output (dBm)
Pmax = 10^((Pmax-30)/10);
PUA = 20  %user power budget (dBm)
PUA = 10^((PUA-30)/10);
PBS = 40  %BS max power consumption (dBm)
PBS = 10^((PBS-30)/10);
tau = 1e-3; %delay
Tf = 200;
Ts = 500;

sigma1 = 1 %noise power is -90dbm and normalize to 1
sigma2 = 1
sigma3 = 1

scale_factor = sqrt(10^12);

N = 64  %BS attenna number
Nrf = 6  %BS RF chain number
Na = 8  %user A attenna number
Nrfa = 2  %user A RF chain number
Nb = 8  %user B attenna number
Nrfb = 2  %user B RF chain number
D1 = Nrfa  %data symbols
D2 = Nrfb
D3 = Nrfa
Dx = 5;
Dy = 50;
n_clu = 4;
n_ray = 5;
path_num = n_clu*n_ray+1;

belta1 = 3;
belta2 = 3;
belta3 = 2.4;
%% channel parameters initialization
position_BS = [0,0,10];
position_userA = [Dx,Dy,1];
position_userB = [-Dx,Dy,1];

AoA1 = atan((1-10)/sqrt(Dx^2+Dy^2));
AoA2 = atan((10-1)/sqrt(Dx^2+Dy^2));
AoA3 = atan((1-1)/2*Dx);
dist1 = norm(position_userA-position_BS,'fro');
dist2 = norm(position_userB-position_BS,'fro');
dist3 = norm(position_userA-position_userB,'fro');

tau1 = tau;
tau2 = tau*(N*Nb)/(N*Na);
tau3 = tau*(Na*Nb)/(N*Na);

tts_tau1 = tau*(Nrf*Nrfa)/(N*Na);
tts_tau2 = tau*(Nrf*Nrfb)/(N*Na);
tts_tau3 = tau*(Nrfa*Nrfb)/(N*Na);

% [nlos_AOA1,nlos_AOD1] = generate_nLOS_angle(n_clu,n_ray);
% [nlos_AOA2,nlos_AOD2] = generate_nLOS_angle(n_clu,n_ray);
% [nlos_AOA3,nlos_AOD3] = generate_nLOS_angle(n_clu,n_ray);
load('nlos_angle1.mat')
load('nlos_angle2.mat')
load('nlos_angle3.mat')

[H1,H1_tau,H1_tau_tts,Ar1,At1] = mmWavechannel(N,Na,nlos_AOA1,nlos_AOD1,dist1,belta1,AoA1,tau1,tts_tau1,scale_factor);
[H2,H2_tau,H2_tau_tts,Ar2,At2] = mmWavechannel(Nb,N,nlos_AOA2,nlos_AOD2,dist2,belta2,AoA2,tau2,tts_tau2,scale_factor);
[H3,H3_tau,H3_tau_tts,Ar3,At3] = mmWavechannel(Nb,Na,nlos_AOA3,nlos_AOD3,dist3,belta3,AoA3,tau3,tts_tau3,scale_factor);
%% long-term annalog matrix initialization

Theta_Fa = 2*pi*rand(Na,Nrfa);
Theta_U1 = 2*pi*rand(N,Nrf);
Theta_U2 = 2*pi*rand(N,Nrf);
Theta_Fb = 2*pi*rand(Nb,Nrfb);

Fa = exp(1j*Theta_Fa);
U1 = exp(1j*Theta_U1);
U2 = exp(1j*Theta_U2);
Fb = exp(1j*Theta_Fb);

f_Theta_U1 = zeros(size(U1));
f_Theta_U2 = zeros(size(U2));
f_Theta_Fa = zeros(size(Fa));
f_Theta_Fb = zeros(size(Fb));

initial_rho = 0.5;
w1 = -initial_rho/(initial_rho+alpha);
w2 = -initial_rho*alpha/(initial_rho+alpha);
w3 = -(1-initial_rho)*alpha/(initial_rho+alpha);

% [OMP_rate1,OMP_Fa1,OMP_U1] = compute_OMP_rate(H1,H1_tau,At1,Ar1,Nrfa,Nrf,PUA,Pmax,sigma1);
% [OMP_rate2,OMP_U2,OMP_Fb2] = compute_OMP_rate(H2,H2_tau,At2,Ar2,Nrf,Nrfb,PBS,Pmax,sigma2);
% [OMP_rate3,OMP_Fa3,OMP_Fb3] = compute_OMP_rate(H3,H3_tau,At3,Ar3,Nrfa,Nrfb,PUA,Pmax,sigma3);
% 
% Fa = w1/(w1+w2)*OMP_Fa1+w2/(w1+w2)*OMP_Fa3;
% U1 = OMP_U1;
% U2 = OMP_U2;
% Fb = w2/(w2+w3)*OMP_Fb2 + w3/(w2+w3)*OMP_Fb3;
% 
% [CM_rate1,CM_Fa1,CM_U1] = compute_CM_rate(H1,H1_tau,Nrfa,Nrf,PUA,Pmax,sigma1);
% [CM_rate2,CM_U2,CM_Fb2] = compute_CM_rate(H2,H2_tau,Nrf,Nrfb,PBS,Pmax,sigma2);
% [CM_rate3,CM_Fa3,CM_Fb3] = compute_CM_rate(H3,H3_tau,Nrfa,Nrfb,PUA,Pmax,sigma3);
% 
% Fa = w1/(w1+w2)*CM_Fa1+w2/(w1+w2)*CM_Fa3;
% U1 = CM_U1;
% U2 = CM_U2;
% Fb = w2/(w2+w3)*CM_Fb2 + w3/(w2+w3)*CM_Fb3;
% 
% Theta_Fa = angle(Fa);
% Theta_U1 = angle(U1);
% Theta_U2 = angle(U2);
% Theta_Fb = angle(Fb);

for t = 1:Tf
%% short-term optimization
t
short_term_delay = 0;
short_term_capacity = 0;
short_term_rate = 0;
OMP_rate = 0;
CM_rate = 0;
parfor i = 1:Ts
%generate channel
[H1,H1_tau,H1_tau_tts,Ar1,At1] = mmWavechannel(N,Na,nlos_AOA1,nlos_AOD1,dist1,belta1,AoA1,tau1,tts_tau1,scale_factor);
[H2,H2_tau,H2_tau_tts,Ar2,At2] = mmWavechannel(Nb,N,nlos_AOA2,nlos_AOD2,dist2,belta2,AoA2,tau2,tts_tau2,scale_factor);
[H3,H3_tau,H3_tau_tts,Ar3,At3] = mmWavechannel(Nb,Na,nlos_AOA3,nlos_AOD3,dist3,belta3,AoA3,tau3,tts_tau3,scale_factor);
%compute channel rate
[proposed_rate1,heuristic_rate1,Wa1] = pcccp(Nrfa,Na,N,Nrf,PUA,Pmax,sigma1,Fa,U1,H1);
[proposed_rate2,heuristic_rate2,V2] = pcccp(Nrf,N,Nb,Nrfb,PBS,Pmax,sigma2,U2,Fb,H2);
[proposed_rate3,heuristic_rate3,Wa3] = pcccp(Nrfa,Na,Nb,Nrfb,PUA,Pmax,sigma3,Fa,Fb,H3);

%compute offloading ratio and delay (proposed)
k1=L/(B1*proposed_rate1);
kE=L/FE;
k2=alpha*L/(B2*proposed_rate2);
kL=L/FL;
k3=alpha*L/(B3*proposed_rate3);
[proposed_rho,delay,condition]=best_partition(k1,kE,k2,kL,k3);
short_term_delay = short_term_delay + compute_total_delay(k1,kE,k2,kL,k3,proposed_rho);

short_term_rate = short_term_rate + w1*proposed_rate1 + w2*proposed_rate2 + w3*proposed_rate3;

short_term_capacity = short_term_capacity + w1*compute_capacity(H1,U1,Fa,sigma1) + w2*compute_capacity(H2,Fb,U2,sigma2) + w3*compute_capacity(H3,Fb,Fa,sigma3);
%compute offloading ratio and delay (heuristic)
% k1h=L/(B1*heuristic_rate1);
% k2h=alpha*L/(B2*heuristic_rate2);
% k3h=alpha*L/(B3*heuristic_rate3);
% [heuristic_rho,delay,condition]=best_partition(k1h,kE,k2h,kL,k3h);
% heuristic_delay(t) = compute_total_delay(k1h,kE,k2h,kL,k3h,heuristic_rho);

% OMP =======================
%     OMP_rate1 = compute_OMP_rate(H1,H1_tau,At1,Ar1,Nrfa,Nrf,PUA,Pmax,sigma1);
%     OMP_rate2 = compute_OMP_rate(H2,H2_tau,At2,Ar2,Nrf,Nrfb,PBS,Pmax,sigma2);
%     OMP_rate3 = compute_OMP_rate(H3,H3_tau,At3,Ar3,Nrfa,Nrfb,PUA,Pmax,sigma3);
%     OMP_rate = OMP_rate + w1*OMP_rate1+w2*OMP_rate2+w3*OMP_rate3;
% 
%     %single_time_scale CM
%     CM_rate1 = compute_CM_rate(H1,H1_tau,Nrfa,Nrf,PUA,Pmax,sigma1);
%     CM_rate2 = compute_CM_rate(H2,H2_tau,Nrf,Nrfb,PBS,Pmax,sigma2);
%     CM_rate3 = compute_CM_rate(H3,H3_tau,Nrfa,Nrfb,PUA,Pmax,sigma3);
%  CM_rate = CM_rate + w1*CM_rate1+w2*CM_rate2+w3*CM_rate3;    
end
proposed_delay(t) = short_term_delay/Ts;
capacity(t) =  short_term_capacity/Ts;
averge_proposed_rate(t) = short_term_rate/Ts;
% averge_OMP_rate = OMP_rate/Ts;
% averge_CM_rate = CM_rate/Ts;
%% long-term optimization
rho_t = t^(-0.6);
gamma_t = t^(-0.9);
varpi = 0.1;

% if t<=1
% w1 = -proposed_rho*L/B1/proposed_rate1^2;
% w2 = -proposed_rho*alpha*L/B2/proposed_rate2^2;
% w3 = -(1-proposed_rho)*alpha*L/B3/proposed_rate3^2;
% w1 = -w1/(w1+w2+w3);
% w2 = -w2/(w1+w2+w3);
% w3 = -w3/(w1+w2+w3);
% end
Wa1 = eye(Nrfa);
V2 = eye(Nrf);
Wa3 = eye(Nrfb); 
g_Theta_U1 = w1*compute_Theta_U1_gradient(Wa1,Fa,H1,U1,sigma1);
g_Theta_U2 = w2*compute_Theta_U2_gradient(V2,U2,H2,Fb,sigma2);
g_Theta_Fa = w1*compute_Theta_U2_gradient(Wa1,Fa,H1,U1,sigma1) + w3*compute_Theta_U2_gradient(Wa3,Fa,H3,Fb,sigma3);
g_Theta_Fb = w2*compute_Theta_U1_gradient(V2,U2,H2,Fb,sigma2) + w3*compute_Theta_U1_gradient(Wa3,Fa,H3,Fb,sigma3);

% proposed_rate1
gradient_norm(t) = norm(g_Theta_U1,'fro')+norm(g_Theta_Fa,'fro');

f_Theta_U1 = (1-rho_t)*f_Theta_U1 + rho_t*g_Theta_U1;
f_Theta_U2 = (1-rho_t)*f_Theta_U2 + rho_t*g_Theta_U2;
f_Theta_Fa = (1-rho_t)*f_Theta_Fa + rho_t*g_Theta_Fa;
f_Theta_Fb = (1-rho_t)*f_Theta_Fb + rho_t*g_Theta_Fb;

Theta_U1_bar = Theta_U1 - f_Theta_U1/(2*varpi);
Theta_U2_bar = Theta_U2 - f_Theta_U2/(2*varpi);
Theta_Fa_bar = Theta_Fa - f_Theta_Fa/(2*varpi);
Theta_Fb_bar = Theta_Fb - f_Theta_Fb/(2*varpi);

Theta_U1 = (1-gamma_t)*Theta_U1 + gamma_t*Theta_U1_bar;
Theta_U2 = (1-gamma_t)*Theta_U2 + gamma_t*Theta_U2_bar;
Theta_Fa = (1-gamma_t)*Theta_Fa + gamma_t*Theta_Fa_bar;
Theta_Fb = (1-gamma_t)*Theta_Fb + gamma_t*Theta_Fb_bar;

Fa = exp(1j*Theta_Fa);
U1 = exp(1j*Theta_U1);
U2 = exp(1j*Theta_U2);
Fb = exp(1j*Theta_Fb);
end
% plot(proposed_delay,'linewidth',1)
% plot(capacity,'linewidth',1)
% xlabel('Number of iterations')
% ylabel('System delay (ms)')
% set(get(gca,'XLabel'),'Fontsize',12)
% set(get(gca,'YLabel'),'Fontsize',12)
figure()
x_axis = 1:Tf;
[ax,line1,line2]=plotyy(x_axis,-capacity,x_axis,proposed_delay*1e3);
xlabel('number of iterations')
set(line1,'linewidth',1.5);
set(line2,'linewidth',1.5,'linestyle',':');
set(ax(1),'Ylim',[12,22])
set(ax(1),'yTick',[12:2:22])
set(ax(2),'Ylim',[3,5])
set(ax(2),'yTick',[3:0.4:5])

legend([line1,line2],'weighted capacity','delay','fontsize',11);
d1=get(ax(1),'ylabel');
set(d1,'Interpreter','latex','string','Weighted capacity (bps/Hz)','fontsize',12);
d2=get(ax(2),'ylabel');
set(d2,'Interpreter','latex','string','System delay (ms)','fontsize',12);
