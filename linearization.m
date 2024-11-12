%% Equilibrium
box_params = struct();
box_params.m = 5;
box_params.I = 10;
box_params.k_list = [3 3 3 3];
box_params.l0_list = [4 4 4 4];
box_params.P_world = [0 0 5 5 0; 0 5 5 0 0];
box_params.P_box = [2 2 4 4 2; 2 4 4 2 2];
box_params.g = -9.81;
t_in = 0;

my_rate_func = @(V_in) box_rate_func(t_in,V_in,box_params);

V_eq = multi_newton_solver(my_rate_func, [1; 1; 0; 0; 0; 0], true);

%% Linearization

J_approx = approximate_jacobian(my_rate_func, V_eq);

my_linear_rate = @(t_in,V_in) J_approx*(V_in-V_eq);
%my_linear_rate(t_in, V_in);

%% Modal Analysis
%[U_mode, omega_n] = eig(Q)
U_mode = 0;%your code here (use eig)
omega_n = 0;%your code here (use eig)
%small number
epsilon = 0;%your code here
V0 = Veq + epsilon*[Umode;0;0;0];
tspan = 0;%your code here
%run the integration of nonlinear system
% [tlist_nonlinear,Vlist_nonlinear] =...
% your_integrator(my_rate_func,tspan,V0,...);
