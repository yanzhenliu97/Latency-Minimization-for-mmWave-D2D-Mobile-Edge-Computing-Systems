function [AOA,AOD] = generate_nLOS_angle(n_clu,n_ray)
angle_spread = 10; %degree
AOA = zeros(n_clu,n_ray);
AOD = zeros(n_clu,n_ray);
for i = 1:n_clu
      theta_clu_r = pi*rand(1);    
      %theta_clu_r = pi;
      theta_clu_t = pi*rand(1);

    for j =1:n_ray
        AOA(i,j) = theta_clu_r + pi*angle_spread/180*rand(1);
        AOD(i,j) = theta_clu_t + pi*angle_spread/180*rand(1);
    end
end