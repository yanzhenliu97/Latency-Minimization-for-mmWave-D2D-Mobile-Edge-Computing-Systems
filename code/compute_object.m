function [result,penalty1,penalty2,penalty3] = compute_object(H,U,V,Fa,Wa,Z,sigma,rho,Ppa,Vout,P,Pmax)
Na = size(Fa,1);
object1 = trace(Z*((V'*U'*H*Fa*Wa-eye(size(Z)))*(V'*U'*H*Fa*Wa-eye(size(Z)))' + sigma^2*V'*U'*U*V)) -log(det(Z));

penalty1 = 0;
penalty2 = 0;
for i =1:Na
    penalty1 = penalty1 + (Ppa(i)-f(Vout(i),Pmax))^2;
    penalty2 = penalty2 + (norm(Fa(i,:)*Wa,'fro')-Vout(i))^2;
end
penalty3 = (sum(Ppa) - P)^2;   
object2 = 1/(2*rho)*(penalty1 + penalty2 + penalty3);

result = object1 + object2;