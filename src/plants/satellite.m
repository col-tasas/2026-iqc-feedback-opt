%% Satellite dynamics (Hill-Clohessy-Wiltshire)

%% Parameters
mc    = 1000;
a     = 6.6e6;
mu_gv = 3.986e14;
omega = sqrt(mu_gv/a^3);

%% State-space matrices
A21 = [3*omega^2, 0,        0; 
       0,         0,        0; 
       0,         0, -omega^2];

A22 = [0,       2*omega, 0; 
      -2*omega, 0,       0; 
       0,       0,       0];

A = [zeros(3), eye(3);
     A21,        A22];

B  = [zeros(3); 
      1/mc*eye(3)];

C  = [eye(3), zeros(3)];

D  = zeros(3);

Bw = [zeros(3,1); 
      [0; 0.01; 0]];

Dw = zeros(3,1);

%% Dimensions

nx = size(A,1);
nu = size(B,2);
ny = size(C,1);
nw = size(Bw,2);

%% Steady-state mappings via null-space
N    = null([A B]);
Pi_y = [C D] * N;
Pi_u = [zeros(nu,nx) eye(nu)] * N;

%% Initial condition

x0 = [50; -140; 270; 0; 0; 0];
