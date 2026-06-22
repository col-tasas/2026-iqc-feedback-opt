%% System dynamics

A = [-1, -4, -1,  3;
      1, -4, -1, -3;
     -1,  4, -1, -9;
      0,  0,  0, -4];

B  = [0; 1; 0; 1];

C  = [1, -1, 0, -4;
      1,  0,  2,  0];

D  = [0; 0];

Bw = [1; 0; 0; 0];

Dw = [0; 0];

nx = size(A,1);
nu = size(B,2);
ny = size(C,1);
nw = size(Bw,2);

%% Steady-state maps

Pi_xu = -A \ B;
Pi_yu =  D - C*(A\B);
Pi_yw =  Dw - C*(A\Bw);

%% Initial Condition

x0 = 2*rand(nx,1) - 1;
