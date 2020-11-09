function result = compute_Theta_U1_gradient(Wa1,Fa,H1,U1,sigma1)
Y1 = eye(size(U1,1)) + H1*Fa*Wa1*Wa1'*Fa'*H1'*U1*pinv(sigma1^2*U1'*U1)*U1';
dR_conj = 1/sigma1^2*(inv(Y1)-U1*inv(U1'*U1)*U1'*inv(Y1))*H1*Fa*Wa1*Wa1'*Fa'*H1'*U1*inv(U1'*U1);
dR = conj(dR_conj);
result = dR.*(1j*U1) - dR_conj.*(1j*conj(U1));