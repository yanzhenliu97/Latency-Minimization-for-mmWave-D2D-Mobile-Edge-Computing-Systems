function result = compute_capacity(H,U,Fa,sigma)

result = log2(det(eye(size(U,2)) + U'*H*Fa*Fa'*H'*U*pinv(sigma^2*U'*U)));
result = real(result);