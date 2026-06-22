clear;
close all;
rng(0);

%% insert your IQCLAB path here %%
error('First insert your IQClab path')
path_to_iqclab_install = ...
run(path_to_iqclab_install)

addpath('sim')
addpath('plants')

%% Load system

run('plants/academic_example.m')

%% Cost functions: Phi1(u) = 1/2 u'Qu u,  Phi2(y) = 1/2 y'Qy y

w_const = 5;

sqrt_Q = randn(nu);
Qu = sqrt_Q'*sqrt_Q + 1e1*eye(nu);

sqrt_Q = randn(ny);
Qy = sqrt_Q'*sqrt_Q + 1e1*eye(ny);

mu = min(eig(Qu));  Lu = max(eig(Qu));
my = min(eig(Qy));  Ly = max(eig(Qy));

%% Engineering constraints and tuning

E = 1;
F = [0, -1];

nc = size(E,1);

u_min = 0.2;
u_max = 1;

alpha = 0.05;
beta  = 10;
rho   = 1;
leak  = 0.001; % improves numeric stability at the cost of slight optimality

% Select controller mode: 0 = hard constraint, 1 = smooth
smooth_hard = 1;

%% True optimal (via fmincon)

obj = @(u) 0.5*u'*Qu*u + 0.5*(Pi_yu*u + Pi_yw*w_const)'*Qy*(Pi_yu*u + Pi_yw*w_const);
Aeq = E + F*Pi_yu;
beq = -F*Pi_yw*w_const;

u_opt = fmincon(obj, zeros(nu,1), [], [], Aeq, beq, u_min, u_max);
y_opt = Pi_yu*u_opt + Pi_yw*w_const;

%% Generalized plant P

N1 = E + F*Pi_yu;
N2 = beta*(E + F*D);
N3 = alpha*(mu*eye(nu) + my*Pi_yu'*D);
N4 = eye(nu) - N3;

% State: [x; eta; lambda]
A_tilde = [A,                    zeros(nx,nu),  zeros(nx,nc);
           alpha*my*Pi_yu'*C,    zeros(nu,nu),  alpha*N1';
           beta*F*C,             zeros(nc,nu),  -leak*ones(nc,nc)];

Bu  = [B;  N3;  N2];
B_w = [Bw; alpha*my*Pi_yu'*Dw; beta*F*Dw];

Cy    = [zeros(nu,nx), eye(nu), zeros(nu,nc)];
Dypwu = zeros(nu, nu+ny+nu+nw+nu);

% Performance: z = [eta; rho*u]
Cz    = [zeros(nu,nx), eye(nu), zeros(nu,nc);
         zeros(nu, nx+nu+nc)               ];
Dzpwu = [zeros(nu, nu+ny+nu+nw+nu);
         zeros(nu, nu+ny+nu), zeros(nu,nw), rho*eye(nu)];

switch smooth_hard
    case 0  % hard constraint controller
        Bp  = [zeros(nx,nu),    zeros(nx,ny),   -alpha*B;
               alpha*eye(nu),   alpha*Pi_yu',    alpha*N4;
               zeros(nc,nu),    zeros(nc,ny),   -alpha*N2];

        Cq  = [zeros(nu,nx),  zeros(nu,nu),  zeros(nu,nc);
               C,             zeros(ny,nu),  zeros(ny,nc);
               zeros(nu,nx),  zeros(nu,nu),  zeros(nu,nc)];
        Dqp = [zeros(nu,nu),  zeros(nu,ny),  -alpha*eye(nu);
               zeros(ny,nu),  zeros(ny,ny),  -alpha*D;
               zeros(nu,nu),  zeros(nu,ny),  -alpha*eye(nu)];
        Dqu = [eye(nu); D; eye(nu)];
        Dqw = [zeros(nu,nw); Dw; zeros(nu,nw)];

    case 1  % smooth controller
        Bp  = [zeros(nx,nu),   zeros(nx,ny),  zeros(nx,nu);
               alpha*eye(nu),  alpha*Pi_yu',  alpha*eye(nu);
               zeros(nc,nu),   zeros(nc,ny),  zeros(nc,nu)];

        Cq  = [zeros(nu,nx),        zeros(nu,nu),  zeros(nu,nc);
               C,                   zeros(ny,nu),  zeros(ny,nc);
               -alpha*my*Pi_yu'*C,  zeros(nu,nu),  -alpha*N1'];
        Dqp = [zeros(nu, nu+ny+nu);
               zeros(ny, nu+ny+nu);
               -alpha*eye(nu), -alpha*Pi_yu', -alpha*eye(nu)];
        Dqu = [eye(nu); D; N4];
        Dqw = [zeros(nu,nw); Dw; -alpha*my*Pi_yu'*Dw];
end

B_tilde = [Bp, B_w, Bu];
C_tilde = [Cq; Cz; Cy];
D_tilde = [Dqp, Dqw, Dqu;
           Dzpwu;
           Dypwu];

P_tilde = minreal(ss(A_tilde, B_tilde, C_tilde, D_tilde));

np    = size(Bp,2);
nq    = size(Cq,1);
nz    = size(Cz,1);
nmeas = size(Cy,1);
ncont = size(Bu,2);

%% True Delta and reduced plant

Delta_Phi1_true = Qu - mu*eye(nu);
Delta_Phi2_true = Qy - my*eye(ny);

% eliminate trivial nabla Phi_1- and nabla Phi_2-channels
P_tilde_ = lft(blkdiag(Delta_Phi1_true, Delta_Phi2_true), P_tilde, nu+ny, nu+ny);

%% IQC synthesis

normal_cone = iqcdelta('cone', ...
                        InputChannel  = 1:nu, ...
                        OutputChannel = 1:nu, ...
                        LinNonlin     = 'NL', ...
                        SectorBounds  = [0, inf], ...
                        SlopeBounds   = [0, inf]);
Delta_cone  = iqcassign(normal_cone, 'usbsr', ...
                        Length            = 2, ...
                        BasisFunctionType = 1, ...
                        PoleLocation      = -1, ...
                        Odd               = 'no');

options           = struct();
options.subopt    = 1.05;
options.condnr    = 1.01;
options.constants = 1e-6*ones(1,3);
options.Pi11pos   = 1e-5;
options.FeasbRad  = 1e8;

[K_syn, gammas] = fRobsyn(P_tilde_, Delta_cone, [nq-nu-ny, nz, nmeas], [np-nu-ny, nw, ncont], options);

[~, idx] = min(gammas);
K_robsyn = K_syn{idx};

%% Build resolvent filter R(s)

s = tf('s');
R = minreal(-K_robsyn * inv(s*eye(nu) - K_robsyn));

%% Simulate

T_sim = 18;

simout = sim('sim/sim_model_constrained.slx');

y = simout.y;
u = simout.u;

figure()
subplot(2,1,1)
plot(u, 'LineWidth', 1.5)
hold on;  grid on
yline(u_opt, '--')
ylabel('u')
xlabel('t')
title('Input')

subplot(2,1,2)
plot(y, 'LineWidth', 1.5)
hold on;  grid on
yline(y_opt, '--')
ylabel('y')
xlabel('t')
title('Output')

sgtitle('Closed-loop trajectories')
