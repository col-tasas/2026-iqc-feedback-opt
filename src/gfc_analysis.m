clear;
close all;
rng(0);

run(../IQClab/IQClab_install.m)

addpath('sim')
addpath('plants')

%% Load system

run('plants/academic_example.m')

%% Cost functions: Phi1(u) = 1/2 u'Qu u,  Phi2(y) = 1/2 y'Qy y

w_const = 10;

sqrt_Q = randn(nu);
Qu = sqrt_Q'*sqrt_Q + 1e-1*eye(nu);

sqrt_Q = randn(ny);
Qy = sqrt_Q'*sqrt_Q + 1e-1*eye(ny);

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
       zeros(nu,nw)];

Cq  = [zeros(nu,nx),  zeros(nu,nu);
       C,              zeros(ny,nu)];

Cy  = [zeros(nu,nx), eye(nu)];

Cz  = [zeros(nu,nx), eye(nu)];    

Dqp   = zeros(nu+ny, nu+ny);
Dqw   = [zeros(nu,nw); Dw];
Dqu   = [eye(nu); D];

Dzpwu = zeros(nu, nu+ny+nw+nu);        

Dyp   = zeros(nu, nu+ny);
Dyw   = zeros(nu, nw);
Dyu   = zeros(nu);

B_tilde = [Bp, B_w, Bu];
C_tilde = [Cq; Cz; Cy];
D_tilde = [Dqp, Dqw, Dqu;
           Dzpwu;
           Dyp, Dyw, Dyu];

P_tilde = ss(A_tilde, B_tilde, C_tilde, D_tilde);

nz = size(Cz,1);

%% True Delta

Delta_Phi1_true = Qu - mu*eye(nu);
Delta_Phi2_true = Qy - my*eye(ny);
Delta_true      = blkdiag(Delta_Phi1_true, Delta_Phi2_true);

%% IQC analysis: sweep over timescale-separation gain k

k_low  = 1e-3;
k_high = 0.2;

nabla_Phi1 = iqcdelta('nabla_Phi1', ...
                       InputChannel  = 1:nu, ...
                       OutputChannel = 1:nu, ...
                       LinNonlin     = 'NL', ...
                       SectorBounds  = [0, Lu-mu], ...
                       SlopeBounds   = [0, Lu-mu]);
nabla_Phi2 = iqcdelta('nabla_Phi2', ...
                       InputChannel  = nu+1:nu+ny, ...
                       OutputChannel = nu+1:nu+ny, ...
                       LinNonlin     = 'NL', ...
                       SectorBounds  = [0, Ly-my], ...
                       SlopeBounds   = [0, Ly-my]);
Delta_Phi1 = iqcassign(nabla_Phi1, 'usbsr', ...
                        Length            = 1, ...
                        PoleLocation      = -1, ...
                        Odd               = 'no');
Delta_Phi2 = iqcassign(nabla_Phi2, 'usbsr', ...
                        Length            = 2, ...
                        BasisFunctionType = 1, ...
                        PoleLocation      = -1, ...
                        Odd               = 'no');
l2_gain = iqcdelta('perf', ...
                    ChannelClass  = 'P', ...
                    InputChannel  = nu+ny+1:nu+ny+nw, ...
                    OutputChannel = nu+ny+1:nu+ny+nz, ...
                    PerfMetric    = 'L2');

options        = struct();
options.Parser = 'LMIlab';

k_vals     = [];
gamma_vals = [];

for k = linspace(k_low, k_high, 20)
    K    = ss([], [], [], -k*eye(nu));
    M    = lft(P_tilde, K);
    prob = iqcanalysis(M, {Delta_Phi1, Delta_Phi2, l2_gain}, options);

    if ~strcmp(prob.gamma, 'infeasible') && ~all(prob.gamma == -1)
        k_vals     = [k_vals;     k         ];
        gamma_vals = [gamma_vals; prob.gamma ];
    end
end

figure()
semilogy(k_vals, gamma_vals, 'LineWidth', 1.5);
grid on
xlabel('Gain k')
ylabel('l2-gain w \rightarrow z')
title('IQC bound vs timescale gain')

[gamma_opt, idx] = min(gamma_vals);
k_opt = k_vals(idx);
K_opt = ss([], [], [], -k_opt*eye(nu));

gamma_iqc  = gamma_opt;
gamma_true = hinfnorm(lft(Delta_true, lft(P_tilde, K_opt)));

%% Simulate

K     = K_opt;
T_sim = 10;

simout = sim('sim/sim_model.slx');

y = simout.y;
u = simout.u;

figure()
subplot(2,1,1)
plot(u, 'LineWidth', 1.5)
hold on;  grid on
yline(u_opt, '--')
ylabel('u')
xlabel('t')
legend('u', 'u_{opt}')
title('Input')

subplot(2,1,2)
plot(y, 'LineWidth', 1.5)
hold on;  grid on
yline(y_opt, '--')
ylabel('y')
xlabel('t')
legend('y_1', 'y_2', 'y_{opt}')
title('Output')

sgtitle('Closed-loop trajectories')
