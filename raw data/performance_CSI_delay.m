load long_term_convergence_10_14.mat U1 U2 Fa Fb

FE = 1600e6  %edge computing rate
FL = 200e6  %local computing rate
L = 1e6  %bits in a time slot
alpha = 0.02  %compression ratio

B1 = 100e6  %bandwidth
B2 = 100e6
B3 = 100e6
Pmax = 30  %PA max power output (dBm)
Pmax = 10^((Pmax-30)/10);
PUA = 20  %user max power consumption (dBm)
PUA = 10^((PUA-30)/10);
PUA = 50*1e-3;
PBS = 40  %BS max power consumption (dBm)
PBS = 10^((PBS-30)/10);
tau = 1e-3; %delay
Tf = 200;
Ts = 100;

sigma1 = 1 %noise power is -80dbm and normalize to 1
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
%% channel initialization
position_BS = [0,0,10];
position_userA = [Dx,Dy,1];
position_userB = [-Dx,Dy,1];

AoA1 = atan((1-10)/sqrt(Dx^2+Dy^2));
AoA2 = atan((10-1)/sqrt(Dx^2+Dy^2));
AoA3 = atan((1-1)/2*Dx);
dist1 = norm(position_userA-position_BS,'fro');
dist2 = norm(position_userB-position_BS,'fro');
dist3 = norm(position_userA-position_userB,'fro');

% [nlos_AOA1,nlos_AOD1] = generate_nLOS_angle(n_clu,n_ray);
% [nlos_AOA2,nlos_AOD2] = generate_nLOS_angle(n_clu,n_ray);
% [nlos_AOA3,nlos_AOD3] = generate_nLOS_angle(n_clu,n_ray);
load('nlos_angle1.mat')
load('nlos_angle2.mat')
load('nlos_angle3.mat')
%% investigate the effect of the CSI delay 
test_num = 5;
proposed_delay = zeros(test_num,6);
heuristic_delay = zeros(test_num,6);
local_computing_delay = zeros(test_num,6);
edge_computing_delay = zeros(test_num,6);
OMP_delay = zeros(test_num,6);
AO_delay = zeros(test_num,6);
CM_delay = zeros(test_num,6);

for j = 1:6
    tau = (j-2)*1e-3;
for i = 1:test_num
    tau1 = tau;
    tau2 = tau*(N*Nb)/(N*Na);
    tau3 = tau*(Na*Nb)/(N*Na);
    
    tts_tau1 = tau*(Nrf*Nrfa)/(N*Na);
    tts_tau2 = tau*(Nrf*Nrfb)/(N*Na);
    tts_tau3 = tau*(Nrfa*Nrfb)/(N*Na);
    
[H1,H1_tau,H1_tau_tts,Ar1,At1] = mmWavechannel(N,Na,nlos_AOA1,nlos_AOD1,dist1,belta1,AoA1,tau1,tts_tau1,scale_factor);
[H2,H2_tau,H2_tau_tts,Ar2,At2] = mmWavechannel(Nb,N,nlos_AOA2,nlos_AOD2,dist2,belta2,AoA2,tau2,tts_tau2,scale_factor);
[H3,H3_tau,H3_tau_tts,Ar3,At3] = mmWavechannel(Nb,Na,nlos_AOA3,nlos_AOD3,dist3,belta3,AoA3,tau3,tts_tau3,scale_factor);
    %compute channel rate
    [proposed_rate1,heuristic_rate1,Wa1] = pcccp(Nrfa,Na,N,Nrf,PUA,Pmax,sigma1,Fa,U1,H1_tau_tts);
    [proposed_rate2,heuristic_rate2,V2] = pcccp(Nrf,N,Nb,Nrfb,PBS,Pmax,sigma2,U2,Fb,H2_tau_tts);
    [proposed_rate3,heuristic_rate3,Wa3] = pcccp(Nrfa,Na,Nb,Nrfb,PUA,Pmax,sigma3,Fa,Fb,H3_tau_tts);
    proposed_rate1_tau = compute_rate(H1,U1,Fa,Wa1,sigma1);
    proposed_rate2_tau = compute_rate(H2,Fb,U2,V2,sigma2);
    proposed_rate3_tau = compute_rate(H3,Fb,Fa,Wa3,sigma3);
    
    %compute offloading ratio and delay (proposed)
    k1=L/(B1*proposed_rate1_tau);
    kE=L/FE;
    k2=alpha*L/(B2*proposed_rate2_tau);
    kL=L/FL;
    k3=alpha*L/(B3*proposed_rate3_tau);
    [proposed_rho,delay,condition]=best_partition(k1,kE,k2,kL,k3);
    proposed_delay(i,j) =  compute_total_delay(k1,kE,k2,kL,k3,proposed_rho);
    
    %compute offloading ratio and delay (heuristic)
    [Wa1,heuristic_rate1] = heuristic(H1_tau_tts,Fa,U1,sigma1,PUA,Nrfa,Pmax);
    [Wa2,heuristic_rate2] = heuristic(H2_tau_tts,U2,Fb,sigma2,PBS,Nrfb,Pmax);
    [Wa3,heuristic_rate3] = heuristic(H3_tau_tts,Fa,Fb,sigma3,PUA,min(Nrfa,Nrfb),Pmax);
    
    heuristic_rate1_tau = compute_rate(H1,U1,Fa,Wa1,sigma1);
    heuristic_rate2_tau = compute_rate(H2,Fb,U2,V2,sigma2);
    heuristic_rate3_tau = compute_rate(H3,Fb,Fa,Wa3,sigma3);
    
    k1h=L/(B1*heuristic_rate1_tau);
    k2h=alpha*L/(B2*heuristic_rate2_tau);
    k3h=alpha*L/(B3*heuristic_rate3_tau);
    [heuristic_rho,delay,condition]=best_partition(k1h,kE,k2h,kL,k3h);
    heuristic_delay(i,j) = compute_total_delay(k1h,kE,k2h,kL,k3h,heuristic_rho);
    
    %compute offloading ratio and delay (local computing)
    local_computing_delay(i,j) = compute_total_delay(k1,kE,k2,kL,k3,0);
    
    %compute offloading ratio and delay (edge computing)
    edge_computing_delay(i,j) = compute_total_delay(k1,kE,k2,kL,k3,1);
    
    %single_time_scale OMP
    OMP_rate1 = compute_OMP_rate(H1,H1_tau,At1,Ar1,Nrfa,Nrf,PUA,Pmax,sigma1);
    OMP_rate2 = compute_OMP_rate(H2,H2_tau,At2,Ar2,Nrf,Nrfb,PBS,Pmax,sigma2);
    OMP_rate3 = compute_OMP_rate(H3,H3_tau,At3,Ar3,Nrfa,Nrfb,PUA,Pmax,sigma3);
    
    k1O=L/(B1*OMP_rate1);
    k2O=alpha*L/(B2*OMP_rate2);
    k3O=alpha*L/(B3*OMP_rate3);
    [OMP_rho,delay,condition]=best_partition(k1O,kE,k2O,kL,k3O);
    OMP_delay(i,j) = compute_total_delay(k1O,kE,k2O,kL,k3O,OMP_rho);
    
    %single_time_scale AO
%     AO_rate1 = compute_AO_rate(H1,H1_tau,Nrfa,Nrf,PUA,Pmax,sigma1);
%     AO_rate2 = compute_AO_rate(H2,H2_tau,Nrf,Nrfb,PBS,Pmax,sigma2);
%     AO_rate3 = compute_AO_rate(H3,H3_tau,Nrfa,Nrfb,PUA,Pmax,sigma3);
%     
%     k1A=L/(B1*AO_rate1);
%     k2A=alpha*L/(B2*AO_rate2);
%     k3A=alpha*L/(B3*AO_rate3);
%     [AO_rho,delay,condition]=best_partition(k1A,kE,k2A,kL,k3A);
%     AO_delay(i,j) = compute_total_delay(k1A,kE,k2A,kL,k3A,AO_rho);
    
    %single_time_scale CM
    CM_rate1 = compute_CM_rate(H1,H1_tau,Nrfa,Nrf,PUA,Pmax,sigma1);
    CM_rate2 = compute_CM_rate(H2,H2_tau,Nrf,Nrfb,PBS,Pmax,sigma2);
    CM_rate3 = compute_CM_rate(H3,H3_tau,Nrfa,Nrfb,PUA,Pmax,sigma3);
    
    k1C=L/(B1*CM_rate1);
    k2C=alpha*L/(B2*CM_rate2);
    k3C=alpha*L/(B3*CM_rate3);
    [CM_rho,delay,condition]=best_partition(k1C,kE,k2C,kL,k3C);
    CM_delay(i,j) = compute_total_delay(k1C,kE,k2C,kL,k3C,CM_rho);
    
end
end
plot(mean(proposed_delay,1))
hold on
plot(mean(heuristic_delay,1))
hold on 
plot(min(mean(local_computing_delay,1),mean(edge_computing_delay,1)))
hold on 
hold on
plot(mean(AO_delay,1))
hold on 
plot(mean(OMP_delay,1))
hold on 
plot(mean(CM_delay,1))
legend('proposed','heuristic','binary','AO','OMP','CM')