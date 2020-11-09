function [W,rate] = waterfilling(H,power,sigma,D)
[Nr,Nt] = size(H);
[U,lambda,V] = svd(H);
chpower = diag(lambda).^2;
num = sum(chpower>0);
for i =1:num
    chpower(i) = 1/chpower(i);
end

water_bound = zeros(num-1,1);
water_power = zeros(Nt,1);

for i =1:num-1
    water_bound(i) = chpower(i+1)*i-sum(chpower(1:i));
end
water_index = sum(water_bound<power)+1;


water_level = (sum(chpower(1:water_index))+power)/water_index;

for i =1:water_index
    water_power(i) = water_level - chpower(i);
end
%sum(water_power)
%beamforming matrix
W = V(:,1:D)*sqrt(diag(water_power(1:D)));
rate = log(det(eye(Nr) + (1/sigma^2)*H*W*W'*H'));