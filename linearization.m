%% Equilibrium
clear;
box_params = get_box_params();
t_in = 0;
V_in = [1; 1; 0; 0; 0; 0];

my_rate_func = @(V_in) box_rate_func(t_in,V_in,box_params);

V_eq = multi_newton_solver(my_rate_func, V_in, true);

%% Linearization

J_approx = approximate_jacobian(my_rate_func, V_eq);
my_rate_func = @(t_in,V_in) box_rate_func(t_in,V_in,box_params);

my_linear_rate = @(t_in,V_in) J_approx*(V_in-V_eq);
my_linear_rate(t_in, V_in);

DormandPrince = make_DP_struct();
h_ref = 0.1;
p=3;
error_desired=0.0001;
tspan = [0,10];
epsilon = 0.05; %small number

V0 = V_eq + epsilon*box_rate_func(t_in,V_in,box_params);

[tlist_nonlinear,Vlist_nonlinear,~,~,~] =...
explicit_RK_variable_step_integration(my_rate_func,tspan,V0,h_ref,DormandPrince,p,error_desired);
[tlist_linear,Vlist_linear,~,~,~] =...
explicit_RK_variable_step_integration(my_linear_rate,tspan,V0,h_ref,DormandPrince,p,error_desired);

%% PLOT
% figure(1);
% plot(Vlist_nonlinear(1,:),Vlist_nonlinear(2,:),"b-");
% hold on;
% plot(Vlist_linear(1,:),Vlist_linear(2,:),"r--");

%% Modal Analysis
Q = -J_approx(4:6,1:3);
[U_mode, omega_n] = eig(Q);

V0 = V_eq + epsilon*[U_mode;0;0;0];

%run the integration of nonlinear system
[tlist_nonlinear,Vlist_nonlinear,~,~,~] =...
explicit_RK_variable_step_integration(my_rate_func,tspan,V0,h_ref,DormandPrince,p,error_desired);
