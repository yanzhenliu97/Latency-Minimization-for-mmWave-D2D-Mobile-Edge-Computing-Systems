function result = compute_rate(H,U,Fa,Wa,sigma)

result = log2(det(eye(size(U,2)) + U'*H*Fa*Wa*Wa'*Fa'*H'*U*pinv(sigma^2*U'*U)));
result = real(result);


