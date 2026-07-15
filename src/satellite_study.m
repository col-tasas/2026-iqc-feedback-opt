clear;
close all;
rng(0);

run('../IQClab/IQClab_install.m') % <-- make sure to set all paths correctly inside here

addpath('sim')
addpath('plants')

%% Plants

run('plants/satellite.m')

%%  Cost function

w_const = ones(nw,1);
y_ref   = [100; 100; 100];

Qu = 0.1*eye(nu);
Qy = 1*eye(ny);

mu = min(eig(Qu));  Lu = max(eig(Qu));
my = min(eig(Qy));  Ly = max(eig(Qy));

%% Tuning

alpha = 0.0001;
beta  = 0;
rho   = 1;
leak  = 0;


smooth_hard = 0;

% NOTE: Setting smooth_hard to 1 requires to uncomment "R (hard)" in
%       Simulink (and to comment "R (smooth)" to avoid internal instability.

%% Optimal solution via null-space parameterization

xup    = [A B] \ (-Bw * w_const);
u_part = xup(nx+1:end);
y_part = [C D] * xup + Dw * w_const;

obj_alpha = @(a) 0.5*(Pi_u*a + u_part)'*Qu*(Pi_u*a + u_part) ...
                + 0.5*(Pi_y*a + y_part - y_ref)'*Qy*(Pi_y*a + y_part - y_ref);
a_opt = fmincon(obj_alpha, zeros(size(N,2),1), [], [], [], [], [], []);

u_opt  = Pi_u * a_opt;
y_opt  = Pi_y * a_opt;
u_opt2 = Pi_u * a_opt + u_part;
y_opt2 = Pi_y * a_opt + y_part;

%% Engineering constraints (empty for satellite)

E  = zeros(1,nu);
F  = zeros(1,ny);
nc = size(E,1);

%% Generalized plant P

N1 = E + F*Pi_y;
N2 = beta*(E + F*D);
N3 = alpha*(mu*Pi_u' + my*Pi_y'*D);
N4 = eye(nu) - N3;

A_tilde = [A,                zeros(nx,nu),  zeros(nx,nc);
           alpha*my*Pi_y'*C, zeros(nu,nu),  alpha*N1';
           beta*F*C,         zeros(nc,nu),  zeros(nc,nc)];

Bu  = [B;  N3;  N2];
B_w = [Bw; alpha*my*Pi_y'*Dw; beta*F*Dw];

Cy    = [zeros(nu,nx), eye(nu), zeros(nu,nc)];
Dypwu = zeros(nu, nu+ny+nu+nw+nu);

switch smooth_hard
    case 0  % hard constraint controller
        Bp  = [zeros(nx,nu),   zeros(nx,ny),  -alpha*B;
               alpha*Pi_u',    alpha*Pi_y',    alpha*N4;
               zeros(nc,nu),   zeros(nc,ny),  -alpha*N2];

        Cq  = [zeros(nu,nx),  zeros(nu,nu),  zeros(nu,nc);
               C,             zeros(ny,nu),  zeros(ny,nc);
               zeros(nu,nx),  zeros(nu,nu),  zeros(nu,nc)];
        Dqp = [zeros(nu,nu),  zeros(nu,ny),  -alpha*eye(nu);
               zeros(ny,nu),  zeros(ny,ny),  -alpha*D;
               zeros(nu,nu),  zeros(nu,ny),  -alpha*eye(nu)];
        Dqu = [eye(nu); D; eye(nu)];
        Dqw = [zeros(nu,nw); Dw; zeros(nu,nw)];

    case 1  % continuous controller
        Bp  = [zeros(nx,nu),   zeros(nx,ny),  zeros(nx,nu);
               alpha*Pi_u',    alpha*Pi_y',   alpha*eye(nu);
               zeros(nc,nu),   zeros(nc,ny),  zeros(nc,nu)];

        Cq  = [zeros(nu,nx),       zeros(nu,nu),  zeros(nu,nc);
               C,                  zeros(ny,nu),  zeros(ny,nc);
               -alpha*my*Pi_y'*C,  zeros(nu,nu),  -alpha*N1'];
        Dqp = [zeros(nu, nu+ny+nu);
               zeros(ny, nu+ny+nu);
               -alpha*Pi_u', -alpha*Pi_y', -alpha*eye(nu)];
        Dqu = [eye(nu); D; N4];
        Dqw = [zeros(nu,nw); Dw; -alpha*my*Pi_y'*Dw];
end


% z = [y; eta; rho*u]
Cz    = [C, zeros(ny,nu), zeros(ny,nc);
         zeros(nu,nx), 10*eye(nu), zeros(nu,nc);
         zeros(nu, nx+nu+nc)];
Dzpwu = [zeros(ny, nu+ny+nu), Dw, D;
         zeros(nu, nu+ny+nu+nw+nu);
         zeros(nu, nu+ny+nu+nw), rho*eye(nu)];

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
options.constants = 1e-6*ones(1,3);
options.Pi11pos   = 1e-5;
options.FeasbRad  = 1e8;
options.Parser    = 'Yalmip';
options.Solver    = 'Mosek';

[K_syn, gammas] = fRobsyn(P_tilde_, Delta_cone, [nq-nu-ny, nz, nmeas], [np-nu-ny, nw, ncont], options);

[~, idx] = min(gammas);
K_robsyn = K_syn{idx};

%% Build resolvent filter R(s)

s = tf('s');
R = minreal(-K_robsyn / (s*eye(nu) - K_robsyn));

%% Simulate


u_min = 1*[-20, -20, -20];
u_max = 1*[ 20,  20,  20];
T_sim = 3000;

simout = sim('sim/sat.slx');

t = simout.u.Time;
u = simout.u.Data;
y = simout.y.Data;

% step disturbance at t=1500: switch setpoint
idx_sp = t >= 1500;
u_sp   = repmat(u_opt(:).', numel(t), 1);
y_sp   = repmat(y_opt(:).', numel(t), 1);
u_sp(idx_sp,:) = repmat(u_opt2(:).', sum(idx_sp), 1);
y_sp(idx_sp,:) = repmat(y_opt2(:).', sum(idx_sp), 1);

figure()
hold on
plot(t, u(:,1), 'LineWidth', 1.5)
plot(t, u(:,2), 'LineWidth', 1.5)
plot(t, u(:,3), 'LineWidth', 1.5)
grid on
plot(t, u_sp, '--', 'Color', [0.5 0.5 0.5])
ylabel('u')
xlabel('t')
legend('u_1', 'u_2', 'u_3', 'u_{opt}')
title('Control input')

figure()
hold on
plot(t, y(:,1), 'LineWidth', 1.5)
plot(t, y(:,2), 'LineWidth', 1.5)
plot(t, y(:,3), 'LineWidth', 1.5)
grid on
plot(t, y_sp, '--', 'Color', [0.5 0.5 0.5])
ylabel('y')
xlabel('t')
legend('y_1', 'y_2', 'y_3', 'y_{opt}')
title('Output')

