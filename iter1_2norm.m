function [K, dt_,dest_,f] = iter1(di,dh,fh,const, num_of_sensors)
%parameters
Kk = const.Kk;
I = const.I;
L = const.L;
L1=const.L1;
Rnk = const.Rnk;
Kt = const.Kt;
Dif = const.Dif;
d0 = const.d0;
ll = size(fh,1);
l1 = size(dh,1);
Kl = I*Kk*L*L1;

% update K
cvx_begin quiet
n = size(dh,2);
variable K(l1,ll);
minimize((vec(K*fh-dh))'*vec(K*fh-dh))
cvx_end

% solve for control action using updated K
cvx_begin quiet
variable f(ll,1);
d = K*f;
% minimize (max(abs(di+d)))
minimize((di+d)'*(di+d))
subject to
f<=ones(ll,1)*577;
f>=-ones(ll,1)*577;

cvx_end
%% output

dt_ = I*Kt*L*L1*f; %the real deformation cause by applying f
dest_ = spsensor(ll,I*Kk*L*L1,dt_,num_of_sensors)+randn*0.0078/3; % the measured deformation.

end