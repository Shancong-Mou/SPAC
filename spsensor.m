function dest = spsensor(ll,Kk,dt, num_sensors)

r = ll;
[U, ~,~] = svd(Kk);
Phir = U(:,1:r);
pt = Phir';
A = pt;
[Q,R,P] = qr(A'*A,'vec');
C = zeros(num_sensors,size(dt,1));
for i =1:num_sensors
    C(i,P(i))=1;
end
% sample
Y = C*dt;
% plus error 
Y = Y+ randn(size(Y))*0.0078/3;
% estimate a
a = pinv(C*Kk)*Y;
% reconstruct y
dest = Kk*a;
end