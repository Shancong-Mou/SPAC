% This is code for case study of "SPAC: Sparse Sensor Placement Based Adaptive Control for High Precision Fuselage Assembly"
% Authors:Shancong Mou, Michael Biehler, Xiaowei Yue, Jeffrey H. Hunt, and Jianjun Shi
% Notice that the Data used in this code are not shared, please use your own data to reproduce the experiments 
% Dummy data is used here and this code does not reproduce the results in the paper

clc;clear;
%% initialization 
InitDef=[]; % initial deformation after applying fixture and gravity forces
InitDef1=[];% initial deformation vefore applying fixture and gravity forces
Deideal = []; % ideal result
Detrue = []; % true result
Dbest = []; % best result
Defall_t_Adpt = []; % adpative control result
Defall_est_Adpt = [];% adaptive estimated control result
Defall_t = []; % sparse sensing control result
Defall_est = [];% sparse sensing estimated control result
num_of_sensors = 20;
num_fuselages = 19; % load 19 fuselages
for ii =1: num_fuselages
%% load fuselage parameters
disp(sprintf('---------- Fuselage # %d ----------',ii))
disp(sprintf('  Load fuselage parameters',ii))
% Load designed stiffness matrix K. 
% Designed K is the same for all fuselages.In our experiment, its dimension 
% is dim_stiffxdim_stiff. For demonstration purpose, we use a 600x600 matrix.
dim_stiff = 600; %  dimension of stiffness matrix
K0 = randn(dim_stiff, dim_stiff); %csvread('StiffMat_full_.csv'); 
K0 = 0.5*(K0+K0'); 
K0 = K0 + dim_stiff*eye(dim_stiff);
% load loadvector: load loadvector for the ii-th fuselage
g = randn(dim_stiff,1); % Loadvec = csvread('g.csv');
% Load true stiffness matrix K_t for the ii-th fuselage 
K_t = randn(dim_stiff, dim_stiff);%csvread('StiffMat_true_.csv');
K_t = 0.5*(K_t+K_t'); 
K_t = K_t + dim_stiff*eye(dim_stiff);
% Load shape deformation before applying fixture and gravity forces
d0  = randn(dim_stiff, 1);% csvread('init.csv');
% record maximum initial deformation before applying fixture and gravity
% forces in InitDef1
InitDef1 = [InitDef1; max(abs(d0))]; 
% Find the boundary points: actuators only place on boundary 
ind = [1:62];
% Load the fixture points: please use your own fixture design information
num_fixture = 8; 
fix = sort([1 2 3 4 5 6 7 8]); 
% delete fixture points from both K and f
fix_expand = sort([6*(fix-1)+1 6*(fix-1)+2 6*(fix-1)+3]);
ls = 1:dim_stiff/6;
ls1 = 1:dim_stiff;
ls = ls';
ls1 = ls1';
ls(fix)=[];
ls1(fix_expand)=[];
% extract boundary point's force
indd = [];%
for i =1:size(ind,2)
    for j =1:3
        indd=[indd; (ind(i)-1)*6+j];
    end
end
[~,idx] = ismember(indd,ls1,'rows');
idx(idx==0)=[];
I = zeros(dim_stiff - 3*num_fixture, dim_stiff - 3*num_fixture);
for i =1:size(idx,1)
I(idx(i),idx(i))=1;
end

%% Calculate initial deformation loading the fuselage onto fixtures
% disp(sprintf('Fuselage # %d: calculate initial deformation loading the fuselage onto fixtures.',ii))

% Extract fixture locating points' initial deviation
uinit = -d0(fix_expand);
% Delete fixture points from both K and f
K1=K0;
K1(fix_expand,:)=[];
K1(:,fix_expand)=[];
K2 = K_t;
K2(fix_expand,:)=[];
K2(:,fix_expand)=[];
g(fix_expand)=[];
d0(fix_expand)=[];
Kk = K1\eye(size(K1));
Kt = K2\eye(size(K2));
% apply fixutre deviation and gravity force and update d0
Kcs=K_t;
Kcs(fix_expand,:)=[];
Kcs(:,setdiff([1:dim_stiff],fix_expand))=[];
d0 = Kt*(-Kcs*uinit+g)+d0;
%%%
I( ~any(I,2), : ) = [];  %rows
%
d0 = I*d0;
%
L = zeros(dim_stiff - 3*num_fixture,size(idx,1));
for i =1:size(idx,1)
    L(idx(i),i)=1;
end
De=[];
%
[~,~,rnk] = unique(ind);
Rnk = zeros(size(ind,1),size(ind,1));
for i =1:size(ind,1)
Rnk(i,rnk(i))=1;
end
% finite difference matrix
Dif = eye(size(ind,1));
for i = 1:size(ind,1)-1
    Dif(i,i+1)=-1;
end
% record maximum initial deformation after applying fixture and gravity
% forces in InitDef
InitDef=[InitDef; max(abs(d0))];

%% searching algorithm for actuator locations
disp(sprintf('  Search for actuator locations.',ii))

ll = 100;
lam = 0.000001;
while ll>20
cvx_begin quiet
% Solve for actuator location
variable f(size(idx,1),1);
d1 = I*Kk*L*f;
minimize (max(abs(d1+d0))+lam*norm(f,1)) 
subject to
f <= 577*ones(size(idx,1),1);
f >= -577*ones(size(idx,1),1);
cvx_end
%
f(abs(f)<10e-4)=0;
f0=f;
indxx = find(f~=0);
in1 = find(f==0);
% 
ll = size(indxx,1); % # of actuators
L1 = zeros(size(idx,1),ll);
for i = 1:ll
    L1(indxx(i),i)=1;
end
lam = lam*2;
end
%% solve again for optimal force on those actuators, based on designed
% disp(sprintf('Fuselage # %d: solve again for optimal force on those actuators, based on designed.',ii))

% information Kk
cvx_begin quiet
variable f(size(idx,1),1);
d1 = I*Kk*L*f;
minimize ((d1+d0)'*(d1+d0))
f<=577*ones(size(idx,1),1);
f>=-577*ones(size(idx,1),1);
subject to
f(in1) == 0;
cvx_end
%% Solve for real deformation use the force derived from desinged
% information Kk
% disp(sprintf('Fuselage # %d: solve for real deformation use the force derived from desinged',ii))

% Define constant
const.Kk = Kk;
const.I = I;
const.L = L;
const.L1 = L1;
const.Dif = Dif;
const.Rnk = Rnk;
const.Kt = Kt;
const.d0 = d0;
KK0 = I*Kk*L*L1;
Ktt = I*Kt*L*L1;

% initial value
%
cvx_begin quiet
variable f1(ll,1);
d1 = KK0*f1; % estimated deformation
% minimize (max(abs(d1+d0)))
minimize ((d1+d0)'*(d1+d0))

subject to
f1<=ones(ll,1)*577;
f1>=-ones(ll,1)*577;
cvx_end
d1t = Ktt*f1; % true deformation


% Using the real information Kt to solve for the best control performance 
cvx_begin quiet
variable f1(ll,1);
dt = Ktt*f1;
% minimize (max(abs(dt+d0)))
minimize ((dt+d0)'*(dt+d0))
subject to
f1<=ones(ll,1)*577;
f1>=-ones(ll,1)*577;
cvx_end

Deideal = [ Deideal; max(abs(d1+d0))]; % ideal result
Detrue = [Detrue; max(abs(d1t+d0))]; % true result
Dbest = [Dbest;  max(abs(dt+d0))]; % best result
disp(sprintf('  Maximum deformation after control'))
disp(sprintf('    Ideal one-shot control %d mm.', max(abs(d1+d0))))
disp(sprintf('    True one-shot control %d mm.', max(abs(d1t+d0))))
disp(sprintf('    Best one-shot control %d mm.',  max(abs(dt+d0))))

%% Estimation and control
disp(sprintf('  Estimation and control: Adpative and SPAC',ii))
% disp(sprintf('Run 100 times to simulate the measurement error'))
Def_est_Adpt=[];
Def_t_Adpt = [];
Def_est=[];
Def_t = [];
% repeat 100 times to simualte the randomness in measurement error, here we
% use 2 for demonstartion purpose
for jjj = 1:2
fh = 100*randn(ll,100);
dh = I*Kt*L*L1*fh;
%% adptive control
destAdpt = dh+randn(size(dh))*0.0078/3; % measurement
di = d1t+d0; % initial real deformation after applying spare leaarning control 
%dt_ the real deformation after the control
%dest_ the estimated real deformation after the control
[K_est, dt_Adpt,dest_Adpt] = iter1_2norm(di,destAdpt,fh,const, num_of_sensors); % 

Mdeformt_Adpt =  max(abs(dt_Adpt+di));
Mdeformest_Adpt =  max(abs(dest_Adpt+di));

Def_t_Adpt = [Def_t_Adpt; Mdeformt_Adpt];
Def_est_Adpt = [Def_est_Adpt; Mdeformest_Adpt];

%% sparse sensing control
dest = spsensor(ll,KK0,dh, num_of_sensors);
[K_est, dt_,dest_] = iter1_2norm(di,dest,fh,const, num_of_sensors); % 
Mdeformt =  max(abs(dt_+di));
Mdeformest =  max(abs(dest_+di));
Def_t = [Def_t; Mdeformt];
Def_est = [Def_est; Mdeformest];
end
Defall_t_Adpt = [Defall_t_Adpt Def_t_Adpt]; % all adptive control result
Defall_est_Adpt = [Defall_est_Adpt Def_est_Adpt];% all adptive estimated control result
Defall_t = [Defall_t Def_t]; % all sparse sensing control result
Defall_est = [Defall_est Def_est];% all estimated sparse sensing control result
disp(sprintf('  Median maximum deformation after control:') )
disp(sprintf('    Adptive control: %d mm.', median(Def_t_Adpt)) )
disp(sprintf('    Estimated adptive control: %d mm.', median(Def_est_Adpt)) )
disp(sprintf('    SPAC control: %d mm.', median(Def_t)) )
disp(sprintf('    Estimated SPAC control: %d mm.', median(Def_est)) )
disp(' ');
end
%% Fig 4 Comparison among the initial deviation, expected control result, and true control result
save('AdaptiveAndSparseSensing_2norm.mat')
load('AdaptiveAndSparseSensing_2norm.mat')
figure()
Def_mid = median(Defall_t_Adpt)';
% bar([Dbest, Deideal, Detrue, Def_mid])

hb = bar([InitDef1, Deideal, Detrue],0.9)
hb(1).FaceColor = [0.8500 0.3250 0.0980];
hb(2).FaceColor = [0.9290 0.6940 0.1250];
hb(3).FaceColor =[0 0.4470 0.7410];
legend('Initial deviation','Expected control result uisng one-shot control method', 'True control result using one-shot control method')

% legend('Best achievable control result', 'Predicted control result using one-shot control method','True control result using one-shot control method', 'Median of control result of using feedback control')
xlabel('Serial number')
ylabel('Maximum deviation (mm)')
% ylim([0 0.12])
% saveas(gcf,strcat('Contr_Res_f',num2str(l1(ii)),'f',num2str(jj)','.jpg'), 'jpg')
%
set(gca,'FontSize',16);
% applyhatch_plusC(gcf, '/x.', 'rkb');
set(hb,'LineWidth', 0.5);
%
figure()
Def_mid = median(Defall_t_Adpt)';
% bar([Dbest, Deideal, Detrue, Def_mid])

hb = bar([Deideal, Detrue])
hb(1).FaceColor = [0.9290 0.6940 0.1250];
hb(2).FaceColor =[0 0.4470 0.7410] ;
legend('Expected control result uisng one-shot control method', 'True control result using one-shot control method')
% legend('Best achievable control result', 'Predicted control result using one-shot control method','True control result using one-shot control method', 'Median of control result of using feedback control')
xlabel('Serial number')
ylabel('Maximum deviation (mm)')
% ylim([0 0.12])
% saveas(gcf,strcat('Contr_Res_f',num2str(l1(ii)),'f',num2str(jj)','.jpg'), 'jpg')
set(gca,'FontSize',16);
set(hb,'LineWidth', 0.5);
% applyhatch_plusC(gcf, '.x', 'bk');

%% Fig 7. Comparison among the one-shot control, the adaptive control with 186-dimensional measurement, and the SPAC with 20-dimensional measurement 

Def_mid1 = median(Defall_t)';
Def_mid2 = median(Defall_t_Adpt)';
Def_3sig1 =  std(Defall_t)';
Def_3sig2 =  std(Defall_t_Adpt)';
Errbd1 = median(Defall_est);
model_series = [Detrue Def_mid1 Def_mid2 Errbd1']; 
model_error = [zeros(size(Detrue,1),1) 3*Def_3sig1 3*Def_3sig2 zeros(size(Detrue,1),1)]; 
b = bar(model_series, 'grouped');
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(model_series);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.4660 0.6740 0.1880];
b(3).FaceColor = [0.4940 0.1840 0.5560];
b(4).FaceColor = [0.9290 0.6940 0.1250];

% Plot the errorbars
errorbar(x',model_series,model_error,'k','linestyle','none');
hold off 
legend('Result of one-shot control method',  'Median of SPAC result with 20-dimensional measurement','Median of adaptive control result with 186-dimensional measurement','Predicted mean of SPAC result with 20-dimensional measurement','Location','best')
xlabel('Serial number')
ylabel('Maximum deviation (mm)')
set(gca,'FontSize',16);
set(b,'LineWidth', 1);

%% Fig 5. Comparison of the control performance between one-shot control and adaptive control 
figure()
%
Def_mid = median(Defall_t_Adpt)';
Def_3sig =  std(Defall_t_Adpt)';
ErrbdAdpt = median(Defall_est_Adpt);

model_series = [Detrue, Def_mid, ErrbdAdpt' ]; 
model_error = [zeros(size(Detrue,1),1) 3*Def_3sig zeros(size(Detrue,1),1)]; 
b = bar(model_series, 'grouped');
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(model_series);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.4940 0.1840 0.5560];
% b(3).Facecolor = [0.9290 0.6940 0.1250];
% Plot the errorbars
errorbar(x',model_series,model_error,'k','linestyle','none');
hold off
legend('Control result of one-shot control method', 'Control result of adaptive control method','Predicted control result of adaptive control method')
xlabel('Serial number')
ylabel('Maximum deviation (mm)')
set(gca,'FontSize',16);
set(b,'LineWidth', 1);
% hold on
% plot(ErrbdAdpt,'Color',[1 0 0],'LineWidth', 2,'LineStyle','--')
%% Fig 8. Sensitivity analysis of the SPAC method by varying the number of sensors
% we change number of measureing pints from 5 -> 180 and repeat the whole
% process for each number of measureing pints, this can be time consuming
clc;clear;
res_1 = [];
res_2 = [];
res_3 = [];
disp(sprintf('We change number of measureing pints from 5 -> 180 \n and repeat the whole process for each number of measureing point \n this can be time consuming'))

for num_of_sensors = 5:5:180
disp(sprintf('Number of sensors %d',num_of_sensors))
%% initialization 
InitDef=[]; % initial deformation after applying fixture and gravity forces
InitDef1=[];% initial deformation vefore applying fixture and gravity forces
Deideal = []; % ideal result
Detrue = []; % true result
Dbest = []; % best result
Defall_t_Adpt = []; % adpative control result
Defall_est_Adpt = [];% adaptive estimated control result
Defall_t = []; % sparse sensing control result
Defall_est = [];% sparse sensing estimated control result
num_of_sensors = 20;
num_fuselages = 19; % load 19 fuselages
for ii =1: 19
%% load fuselage parameters
disp(sprintf('---------- Fuselage # %d ----------',ii))
% disp(sprintf('  Load fuselage parameters',ii))
% Load designed stiffness matrix K. 
% Designed K is the same for all fuselages.In our experiment, its dimension 
% is dim_stiffxdim_stiff. For demonstration purpose, we use a 600x600 matrix.
dim_stiff = 600; %  dimension of stiffness matrix
K0 = randn(dim_stiff, dim_stiff); %csvread('StiffMat_full_.csv'); 
K0 = 0.5*(K0+K0'); 
K0 = K0 + dim_stiff*eye(dim_stiff);
% load loadvector: load loadvector for the ii-th fuselage
g = randn(dim_stiff,1); % Loadvec = csvread('g.csv');
% Load true stiffness matrix K_t for the ii-th fuselage 
K_t = randn(dim_stiff, dim_stiff);%csvread('StiffMat_true_.csv');
K_t = 0.5*(K_t+K_t'); 
K_t = K_t + dim_stiff*eye(dim_stiff);
% Load shape deformation before applying fixture and gravity forces
d0  = randn(dim_stiff, 1);% csvread('init.csv');
% record maximum initial deformation before applying fixture and gravity
% forces in InitDef1
InitDef1 = [InitDef1; max(abs(d0))]; 
% Find the boundary points: actuators only place on boundary 
ind = [1:62];
% Load the fixture points: please use your own fixture design information
num_fixture = 8; 
fix = sort([1 2 3 4 5 6 7 8]); 
% delete fixture points from both K and f
fix_expand = sort([6*(fix-1)+1 6*(fix-1)+2 6*(fix-1)+3]);
ls = 1:dim_stiff/6;
ls1 = 1:dim_stiff;
ls = ls';
ls1 = ls1';
ls(fix)=[];
ls1(fix_expand)=[];
% extract boundary point's force
indd = [];%
for i =1:size(ind,2)
    for j =1:3
        indd=[indd; (ind(i)-1)*6+j];
    end
end
[~,idx] = ismember(indd,ls1,'rows');
idx(idx==0)=[];
I = zeros(dim_stiff - 3*num_fixture, dim_stiff - 3*num_fixture);
for i =1:size(idx,1)
I(idx(i),idx(i))=1;
end

%% Calculate initial deformation loading the fuselage onto fixtures
% disp(sprintf('Fuselage # %d: calculate initial deformation loading the fuselage onto fixtures.',ii))

% Extract fixture locating points' initial deviation
uinit = -d0(fix_expand);
% Delete fixture points from both K and f
K1=K0;
K1(fix_expand,:)=[];
K1(:,fix_expand)=[];
K2 = K_t;
K2(fix_expand,:)=[];
K2(:,fix_expand)=[];
g(fix_expand)=[];
d0(fix_expand)=[];
Kk = K1\eye(size(K1));
Kt = K2\eye(size(K2));
% apply fixutre deviation and gravity force and update d0
Kcs=K_t;
Kcs(fix_expand,:)=[];
Kcs(:,setdiff([1:dim_stiff],fix_expand))=[];
d0 = Kt*(-Kcs*uinit+g)+d0;
%%%
I( ~any(I,2), : ) = [];  %rows
%
d0 = I*d0;
%
L = zeros(dim_stiff - 3*num_fixture,size(idx,1));
for i =1:size(idx,1)
    L(idx(i),i)=1;
end
De=[];
%
[~,~,rnk] = unique(ind);
Rnk = zeros(size(ind,1),size(ind,1));
for i =1:size(ind,1)
Rnk(i,rnk(i))=1;
end
% finite difference matrix
Dif = eye(size(ind,1));
for i = 1:size(ind,1)-1
    Dif(i,i+1)=-1;
end
% record maximum initial deformation after applying fixture and gravity
% forces in InitDef
InitDef=[InitDef; max(abs(d0))];

%% searching algorithm for actuator locations
% disp(sprintf('  Search for actuator locations.',ii))

ll = 100;
lam = 0.000001;
while ll>20
cvx_begin quiet
% Solve for actuator location
variable f(size(idx,1),1);
d1 = I*Kk*L*f;
minimize (max(abs(d1+d0))+lam*norm(f,1)) 
subject to
f <= 577*ones(size(idx,1),1);
f >= -577*ones(size(idx,1),1);
cvx_end
%
f(abs(f)<10e-4)=0;
f0=f;
indxx = find(f~=0);
in1 = find(f==0);
% 
ll = size(indxx,1); % # of actuators
L1 = zeros(size(idx,1),ll);
for i = 1:ll
    L1(indxx(i),i)=1;
end
lam = lam*2;
end
%% solve again for optimal force on those actuators, based on designed
% disp(sprintf('Fuselage # %d: solve again for optimal force on those actuators, based on designed.',ii))

% information Kk
cvx_begin quiet
variable f(size(idx,1),1);
d1 = I*Kk*L*f;
minimize ((d1+d0)'*(d1+d0))
f<=577*ones(size(idx,1),1);
f>=-577*ones(size(idx,1),1);
subject to
f(in1) == 0;
cvx_end
%% Solve for real deformation use the force derived from desinged
% information Kk
% disp(sprintf('Fuselage # %d: solve for real deformation use the force derived from desinged',ii))

% Define constant
const.Kk = Kk;
const.I = I;
const.L = L;
const.L1 = L1;
const.Dif = Dif;
const.Rnk = Rnk;
const.Kt = Kt;
const.d0 = d0;
KK0 = I*Kk*L*L1;
Ktt = I*Kt*L*L1;

% initial value
%
cvx_begin quiet
variable f1(ll,1);
d1 = KK0*f1; % estimated deformation
% minimize (max(abs(d1+d0)))
minimize ((d1+d0)'*(d1+d0))

subject to
f1<=ones(ll,1)*577;
f1>=-ones(ll,1)*577;
cvx_end
d1t = Ktt*f1; % true deformation


% Using the real information Kt to solve for the best control performance 
cvx_begin quiet
variable f1(ll,1);
dt = Ktt*f1;
% minimize (max(abs(dt+d0)))
minimize ((dt+d0)'*(dt+d0))
subject to
f1<=ones(ll,1)*577;
f1>=-ones(ll,1)*577;
cvx_end

Deideal = [ Deideal; max(abs(d1+d0))]; % ideal result
Detrue = [Detrue; max(abs(d1t+d0))]; % true result
Dbest = [Dbest;  max(abs(dt+d0))]; % best result
% disp(sprintf('  Maximum deformation after control'))
% disp(sprintf('    Ideal one-shot control %d mm.', max(abs(d1+d0))))
% disp(sprintf('    True one-shot control %d mm.', max(abs(d1t+d0))))
% disp(sprintf('    Best one-shot control %d mm.',  max(abs(dt+d0))))

%% Estimation and control
% disp(sprintf('  Estimation and control: Adpative and SPAC',ii))
% disp(sprintf('Run 100 times to simulate the measurement error'))
Def_est_Adpt=[];
Def_t_Adpt = [];
Def_est=[];
Def_t = [];
% repeat 100 times to simualte the randomness in measurement error, here we
% use 2 for demonstartion purpose
for jjj = 1:2
fh = 100*randn(ll,100);
dh = I*Kt*L*L1*fh;
%% adptive control
destAdpt = dh+randn(size(dh))*0.0078/3; % measurement
di = d1t+d0; % initial real deformation after applying spare leaarning control 
%dt_ the real deformation after the control
%dest_ the estimated real deformation after the control
[K_est, dt_Adpt,dest_Adpt] = iter1_2norm(di,destAdpt,fh,const, num_of_sensors); % 

Mdeformt_Adpt =  max(abs(dt_Adpt+di));
Mdeformest_Adpt =  max(abs(dest_Adpt+di));

Def_t_Adpt = [Def_t_Adpt; Mdeformt_Adpt];
Def_est_Adpt = [Def_est_Adpt; Mdeformest_Adpt];

%% sparse sensing control
dest = spsensor(ll,KK0,dh, num_of_sensors);
[K_est, dt_,dest_] = iter1_2norm(di,dest,fh,const, num_of_sensors); % 
Mdeformt =  max(abs(dt_+di));
Mdeformest =  max(abs(dest_+di));
Def_t = [Def_t; Mdeformt];
Def_est = [Def_est; Mdeformest];
end
Defall_t_Adpt = [Defall_t_Adpt Def_t_Adpt]; % all adptive control result
Defall_est_Adpt = [Defall_est_Adpt Def_est_Adpt];% all adptive estimated control result
Defall_t = [Defall_t Def_t]; % all sparse sensing control result
Defall_est = [Defall_est Def_est];% all estimated sparse sensing control result
% disp(sprintf('  Median maximum deformation after control:') )
% disp(sprintf('    Adptive control: %d mm.', median(Def_t_Adpt)) )
% disp(sprintf('    Estimated adptive control: %d mm.', median(Def_est_Adpt)) )
% disp(sprintf('    SPAC control: %d mm.', median(Def_t)) )
% disp(sprintf('    Estimated SPAC control: %d mm.', median(Def_est)) )
% disp(' ');
end

Def_mid1 = median(Defall_t)';
Def_mid2 = median(Defall_t_Adpt)';
%%
res_1 = [res_1 mean((Detrue-Def_mid2)./Detrue)];
res_2 = [res_2 mean((Detrue-Def_mid1)./Detrue)];
res_3 = [res_3 mean((Def_mid1-Def_mid2)./Def_mid1)];
end
%
% plot F8
% load('res_1.mat')
% load('res_2.mat')
% load('res_3.mat')
%
axix = 20:5:180;
plot(axix, res_1(4:end),'LineWidth', 2)
hold on
plot(axix, res_2(4:end),'LineStyle', '--','LineWidth',2)
legend('APPI of adaptive control result with 186-dimensional measurement',  'APPI of SPAC result with 20-dimensional measurement','Location','best')
xlabel('Number of sensors')
ylabel('APPI')
set(gca,'FontSize',16);




