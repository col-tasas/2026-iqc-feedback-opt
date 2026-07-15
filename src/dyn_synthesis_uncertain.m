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
Qy = sqrt_Q'*sqrt_Q;

% sector constants
mu = min(eig(Qu));  Lu = max(eig(Qu));
my = min(eig(Qy));  Ly = max(eig(Qy));

% Formulate Objective & Equality Constraints in terms of u
obj = @(u) 0.5*u'*Qu*u + 0.5*(Pi_yu*u + Pi_yw*w_const)'*Qy*(Pi_yu*u + Pi_yw*w_const);

% Find true solution
u_opt = fmincon(obj, zeros(nu,1), [], [], [], [], [], []);
y_opt = Pi_yu * u_opt + Pi_yw*w_const;

%% Nominal reduced model

G_real = ss(A, [B Bw], C, [D Dw]);

rspec = reducespec(G_real, "balanced");
G_red = getrom(rspec, Order=2, Method="MatchDC");

figure()
bode(G_real);  
hold on;  
grid on
bode(G_red)
legend('Original','Reduced')

%% Frequency-domain model error weight W(s)

% Increase the scalar coefficient to model larger truncation errors:
% w_inf = 0.07;  % (order 3);
w_inf = 0.27;  % (order 2);
% w_inf = 0.67;  % (order 1);
% w_inf = 2.85;  % (order 0);
s = tf('s');
W = w_inf * s/(s - max(real(eig(A))));

E_err = minreal(G_real - G_red);

figure()
sigma(E_err);  hold on;  grid on
sigma(W)
legend('E(s)', 'W(s)')

%% Generalized plant P

[Ar, B0r, Cr, D0r] = ssdata(G_red);
Br  = B0r(:,1:nu);       
Bwr = B0r(:,nu+1:end);
Dr  = D0r(:,1:nu);       
Dwr = D0r(:,nu+1:end);

W_full = W * eye(ny);
[Aw, Bwf, Cwf, Dwf] = ssdata(W_full);

[nr,  ~] = size(Ar);
[nwf, ~] = size(Aw);
npd = ny;
nqd = nu + nw;

nx_aug = nr + nwf + nu;

A_tilde = [Ar,              zeros(nr,nwf),   zeros(nr,nu);
           zeros(nwf,nr),   Aw,              zeros(nwf,nu);
           my*Pi_yu'*Cr,    my*Pi_yu'*Cwf,   zeros(nu,nu) ];

Bp_new = [zeros(nr,nu),   zeros(nr,ny);
          zeros(nwf,nu),  zeros(nwf,ny);
          eye(nu),         Pi_yu'          ];

Bpd    = [zeros(nr,npd);
          Bwf;
          my*Pi_yu'*Dwf];

Bw_new = [Bwr;
          zeros(nwf,nw);
          my*Pi_yu'*Dwr];

Bu_new = [Br;
          zeros(nwf,nu);
          mu*eye(nu) + my*Pi_yu'*Dr];

B_tilde = [Bp_new, Bpd, Bw_new, Bu_new];

Cq     = [zeros(nu,nx_aug);
          Cr, Cwf, zeros(ny,nu)];
Dqp    = zeros(nu+ny, nu+ny);
Dqpd   = [zeros(nu,npd); Dwf];
Dqw    = [zeros(nu,nw); Dwr];
Dqu    = [eye(nu); Dr];

Cqd    = zeros(nqd, nx_aug);
Dqdp   = zeros(nqd, nu+ny);
Dqdpd  = zeros(nqd, npd);
Dqdw   = [zeros(nu,nw); eye(nw)];
Dqdu   = [eye(nu); zeros(nw,nu)];

rho = 100;
Cz    = [zeros(nu,nr+nwf), eye(nu);
         zeros(nu,nx_aug)          ];
Dzp   = zeros(2*nu, nu+ny);
Dzpd  = zeros(2*nu, npd);
Dzw   = zeros(2*nu, nw);
Dzu   = [zeros(nu); rho*eye(nu)];

Cy   = [zeros(nu,nr+nwf), eye(nu)];
Dyp  = zeros(nu, nu+ny);
Dypd = zeros(nu, npd);
Dyw  = zeros(nu, nw);
Dyu  = zeros(nu);

C_tilde = [Cq; Cqd; Cz; Cy];
D_tilde = [Dqp,  Dqpd,  Dqw,  Dqu;
           Dqdp, Dqdpd, Dqdw, Dqdu;
           Dzp,  Dzpd,  Dzw,  Dzu;
           Dyp,  Dypd,  Dyw,  Dyu];

P_tilde = minreal(ss(A_tilde, B_tilde, C_tilde, D_tilde));

nz    = size(Cz, 1);
nmeas = size(Cy, 1);
ncont = size(Bu_new, 2);

%% True Delta and reduced plant

Delta_Phi1_true = Qu - mu*eye(nu);
Delta_Phi2_true = Qy - my*eye(ny);

% eliminate trivial nabla Phi_1-channel
P_tilde_ = lft(Delta_Phi1_true, P_tilde, nu, nu);

%% IQC synthesis

nabla_Phi2 = iqcdelta('nabla_Phi2', ...
                       InputChannel  = 1:ny, ...
                       OutputChannel = 1:ny, ...
                       LinNonlin     = 'NL', ...
                       SectorBounds  = [0, Ly-my], ...
                       SlopeBounds   = [0, Ly-my]);
Delta_Phi2 = iqcassign(nabla_Phi2, 'usbsr', ...
                        Length            = 4,  ... %Zames-Falb filter length
                        BasisFunctionType = 1,  ... %Zames-Falb parametrization
                        PoleLocation      = -1, ... %Zames-Falb parametrization
                        Odd               = 'no');

model_err = iqcdelta('model_mismatch', ...
                      InputChannel         = ny+1:ny+nqd, ...
                      OutputChannel        = ny+1:ny+npd, ...
                      StaticDynamic        = 'D', ...
                      Structure            = 'FB', ...
                      NormBounds           = 1);

Delta_mm  = iqcassign(model_err, 'ultid', ...
                       Length            = 3, ...
                       BasisFunctionType = 1, ...
                       PoleLocation      = -1);

options           = struct();
options.subopt    = 1.05;
options.constants = 1e-5*ones(1,3);
options.Pi11pos   = 1e-4;
options.FeasbRad  = 1e6;
options.Parser    = 'Yalmip';
options.Solver    = 'Mosek';

[K_syn, gammas] = fRobsyn(P_tilde_, {Delta_Phi2, Delta_mm}, [ny+nqd, nz, nmeas], [ny+npd, nw, ncont], options);

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
legend('y_1', 'y_2', 'y_{opt}')
title('Output')

sgtitle('Closed-loop trajectories')
