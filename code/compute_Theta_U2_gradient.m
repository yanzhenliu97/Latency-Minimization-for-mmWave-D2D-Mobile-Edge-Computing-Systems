function result = compute_Theta_U2_gradient(V2,U2,H2,Fb,sigma2)
Y2 = eye(size(Fb,2)) + 1/sigma2^2*Fb'*H2*U2*V2*V2'*U2'*H2'*Fb*inv(Fb'*Fb);
dR_conj = 1/sigma2^2*H2'*Fb*inv(Fb'*Fb)*inv(Y2)*Fb'*H2*U2*V2*V2';
dR = conj(dR_conj);
result = dR.*(1j*U2) - dR_conj.*(1j*conj(U2));