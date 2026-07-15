clear;
close all;
rng(0);

run('../IQClab/IQClab_install.m') % <-- make sure to set all paths correctly inside here

addpath('sim')
addpath('plants')

%% Load system

run('plants/academic_example.m')

%% Cost functions: Phi1(u) = 1/2 u'Qu u,  Phi2(y) = 1/2 y'Qy y

w_const = 10;

sqrt_Q = randn(nu);
Qu = sqrt_Q'*sqrt_Q + 1e1*eye(nu);

sqrt_Q = randn(ny);
Qy = sqrt_Q'*sqrt_Q + 1e1*eye(ny);

% sector constants
mu = min(eig(Qu));  Lu = max(eig(Qu));
my = min(eig(Qy));  Ly = max(eig(Qy));

% Formulate Objective & Equality Constraints in terms of u
obj = @(u) 0.5*u'*Qu*u + 0.5*(Pi_yu*u + Pi_yw*w_const)'*Qy*(Pi_yu*u + Pi_yw*w_const);

% Find true solution
u_opt = fmincon(obj, zeros(nu,1), [], [], [], [], [], []);
y_opt = Pi_yu * u_opt + Pi_yw*w_const;

%% Generalized plant P

A_tilde = [A,            zeros(nx,nu);
           my*Pi_yu'*C,  zeros(nu,nu)];

Bp  = [zeros(nx,nu),  zeros(nx,ny);
       eye(nu),        Pi_yu'       ];

Bu  = [B;
       mu*eye(nu) + my*Pi_yu'*D];

B_w = [Bw;
       my*Pi_yu'*Dw];

Cq  = [zeros(nu,nx),  zeros(nu,nu);
       C,              zeros(ny,nu)];

Cy  = [zeros(nu,nx), eye(nu)];

% z = [eta; rho*u]
rho   = 100;
Cz    = [zeros(nu,nx), eye(nu);    
         zeros(nu,nx+nu)          ];

Dzpwu = [zeros(nu, nu+ny+nw+nu);
         zeros(nu, nu+ny+nw), rho*eye(nu)];

Dqp   = zeros(nu+ny, nu+ny);
Dqw   = [zeros(nu,nw); Dw];
Dqu   = [eye(nu); D];

Dyp   = zeros(nu, nu+ny);
Dyw   = zeros(nu, nw);
Dyu   = zeros(nu);

B_tilde = [Bp, B_w, Bu];
C_tilde = [Cq; Cz; Cy];
D_tilde = [Dqp, Dqw, Dqu;
           Dzpwu;
           Dyp, Dyw, Dyu];

P_tilde = ss(A_tilde, B_tilde, C_tilde, D_tilde);

nz    = size(Cz, 1);
nmeas = size(Cy, 1);
ncont = size(Bu, 2);

%% True Delta and reduced plant

Delta_Phi1_true = Qu - mu*eye(nu);
Delta_Phi2_true = Qy - my*eye(ny);

% eliminate trivial nabla Phi_1-channel (scalar case: Delta_Phi1_true = 0 for n_u=1)
P_tilde_ = lft(Delta_Phi1_true, P_tilde, nu, nu);

%% IQC synthesis

nabla_Phi2 = iqcdelta('nabla_Phi2', ...
                     InputChannel  = 1:ny, ...
                     OutputChannel = 1:ny, ...
                     LinNonlin     = 'NL', ...
                     SectorBounds  = [0, Ly-my], ...
                     SlopeBounds   = [0, Ly-my]);
Delta_Phi2 = iqcassign(nabla_Phi2, 'usbsr', ...
                       Length            = 2,  ... %Zames-Falb filter length
                       BasisFunctionType = 1,  ... %Zames-Falb parametrization
                       PoleLocation      = -1, ... %Zames-Falb parametrization
                       Odd               = 'no');

options           = struct();
options.subopt    = 1.05;
options.constants = 1e-6*ones(1,3);
options.Pi11pos   = 1e-5;
options.FeasbRad  = 1e6;

[K_syn, gammas] = fRobsyn(P_tilde_, Delta_Phi2, [ny, nz, nu], [ny, nw, nu], options);

[~, idx] = min(gammas);
K_robsyn = K_syn{idx};


%% Simulate

K     = K_robsyn;
T_sim = 10;

simout = sim('sim/sim_model.slx');

y = simout.y;
u = simout.u;

figure()
subplot(2,1,1)
plot(u, 'LineWidth', 1.5)
hold on
grid on
yline(u_opt, '--')
ylabel('u')
xlabel('t')
legend('u', 'u_{opt}')
title('Input')

subplot(2,1,2)
plot(y, 'LineWidth', 1.5)
hold on
grid on
yline(y_opt, '--')
ylabel('y')
xlabel('t')
legend('y_1', 'y_2', 'y_{opt}')
title('Output')

sgtitle('Closed loop trajectories')
