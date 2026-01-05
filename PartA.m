%% Control Engineering Coursework - MATLAB Code
% Nonlinear system definition and automatic linearization

clear; clc;

%% System Parameters
F_val = 1;      % Friction coefficient
M_val = 1;      % Cart mass
g_val = 9.81;   % Gravitational acceleration
L_val = 1;      % Pendulum length

%% Define Symbolic Variables
syms x1 x2 x3 x4 u F M g L real

% State vector: x = [s; s_dot; phi; phi_dot]
%   x1 = s       (cart position)
%   x2 = s_dot   (cart velocity)
%   x3 = phi     (pendulum angle)
%   x4 = phi_dot (pendulum angular velocity)

%% Define Nonlinear State Equations
% From the equations of motion (with mu = d1 = 0):
%   M*s_ddot + F*s_dot = u
%   phi_ddot - (g/L)*sin(phi) + (1/L)*s_ddot*cos(phi) = 0

% Solving for x_dot:
f1 = x2;                                              % x1_dot = x2
f2 = -F/M * x2 + u/M;                                 % x2_dot = (u - F*x2)/M
f3 = x4;                                              % x3_dot = x4
f4 = (g/L)*sin(x3) + (F/(M*L))*x2*cos(x3) - (u/(M*L))*cos(x3);  % x4_dot

% State equations vector
f = [f1; f2; f3; f4];

disp('Nonlinear State Equations f(x,u):');
disp(f);

%% Compute Jacobians (Linearization)
x = [x1; x2; x3; x4];

% A = df/dx (Jacobian with respect to state)
A_sym = jacobian(f, x);
disp('Symbolic A matrix (df/dx):');
disp(A_sym);

% B = df/du (Jacobian with respect to input)
B_sym = jacobian(f, u);
disp('Symbolic B matrix (df/du):');
disp(B_sym);

%% Evaluate at Equilibrium Point x_e = [0; 0; 0; 0], u_e = 0
% Equilibrium 1: pendulum down (stable)
x_e1 = [0; 0; 0; 0];
u_e = 0;

A_at_eq1 = subs(A_sym, [x1, x2, x3, x4, u, F, M, g, L], ...
                       [0, 0, 0, 0, 0, F_val, M_val, g_val, L_val]);
B_at_eq1 = subs(B_sym, [x1, x2, x3, x4, u, F, M, g, L], ...
                       [0, 0, 0, 0, 0, F_val, M_val, g_val, L_val]);

A = double(A_at_eq1);
B = double(B_at_eq1);

disp('=== Linearized at x_e1 = [0; 0; 0; 0] (pendulum down) ===');
disp('A matrix:'); disp(A);
disp('B matrix:'); disp(B);


%% Output Matrices
C = [1 0 0 0; 0 0 1 0];

%% Controllability Analysis
% Reachability matrix: R = [B, A*B, A^2*B, A^3*B]
R = ctrb(A, B);
disp('Reachability Matrix R:');
disp(R);

% Check rank
rank_R = rank(R);
disp(['Rank of R = ', num2str(rank_R)]);

n = size(A, 1);  % Number of states
if rank_R == n
    disp('System (A, B) is CONTROLLABLE');
else
    disp('System (A, B) is NOT controllable');
end

%% Display in standard state space form
xdot = A * x + B * u;
y = C * x;

%% A5 
C_reg = [1 0 0 0];

% s = 0
H1 = [- A B; C_reg 0];
rank(H1)

% s = j w
syms omega real
s = 1j * omega; 

H2 = [s*eye(4) - A,  B;
     C_reg,         0];

rank(H2)

% s = -jw
H3 = [-s*eye(4) - A,  B;
     C_reg,         0];

rank(H3)

%% A7
S = [0 0 0; 0 0 omega; 0 -omega 0];

P = [0, 0, 0;
     1, 0, 0;
     0, 0, 0;
     -1, 0, 0];

size_A = size(A);
size_S = size(S);

A_aug = [A P; zeros(size_S(1), size_A(2)) S];
C_aug = [C_reg zeros(2,3)];

O = obsv(A_aug, C_aug);

rank(obsv(A_aug, C_aug))

% d1 is constant, so its dynamics are just d1_dot = 0
S_d1 = 0;  % scalar

% Only first column of P (how d1 affects x)
P_d1 = [0; 1; 0; -1];

% 5-state augmented system
A_aug_reduced = [A, P_d1; zeros(1,4), S_d1];
C_aug_reduced = [C, zeros(2,1)];

rank(obsv(A_aug_reduced, C_aug_reduced))