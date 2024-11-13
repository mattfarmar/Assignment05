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
epsilon = 0.25; %small number

V0 = V_eq + epsilon*box_rate_func(t_in,V_in,box_params);

[tlist_nonlinear,Vlist_nonlinear,~,~,~] =...
explicit_RK_variable_step_integration(my_rate_func,tspan,V0,h_ref,DormandPrince,p,error_desired);
[tlist_linear,Vlist_linear,~,~,~] =...
explicit_RK_variable_step_integration(my_linear_rate,tspan,V0,h_ref,DormandPrince,p,error_desired);

figure()
subplot(3,1,1)
plot(tlist_nonlinear, Vlist_nonlinear(1,:), 'b')
hold on
plot(tlist_linear, Vlist_linear(1,:), 'r-')
ylabel('X')
title('X vs. Time (Nonlinear/Linear)')
subplot(3,1,2)
plot(tlist_nonlinear, Vlist_nonlinear(2,:), 'b')
hold on
plot(tlist_linear, Vlist_linear(2,:), 'r-')
ylabel('Y')
title('Y vs. Time (Nonlinear/Linear)')
subplot(3,1,3)
plot(tlist_nonlinear, Vlist_nonlinear(3,:), 'b')
hold on
plot(tlist_linear, Vlist_linear(3,:), 'r-')
ylabel('Theta')
title('Theta vs. Time (Nonlinear/Linear)')
xlabel('Time (s)')
legend('Nonlinear', 'Linear')

%% PLOT
% figure(1);
% plot(Vlist_nonlinear(1,:),Vlist_nonlinear(2,:),"b-");
% hold on;
% plot(Vlist_linear(1,:),Vlist_linear(2,:),"r--");

%% Modal Analysis
Q = -J_approx(4:6,1:3);
[U_mode, omega_n] = eig(Q);


% MODE 1
U_mode1 = U_mode(:,1);
omega_n1 = omega_n(1,1);

V0 = V_eq + epsilon*[U_mode1;0;0;0];

%run the integration of nonlinear system
[tlist_nonlinear,Vlist_nonlinear,~,~,~] =...
explicit_RK_variable_step_integration(my_rate_func,tspan,V0,h_ref,DormandPrince,p,error_desired);

figure()
plot(Vlist_nonlinear(1,:), Vlist_nonlinear(2,:), '-b')
hold on
x_modal1 = V_eq(1) + epsilon*U_mode1(1)*cos(omega_n1*tlist_nonlinear);
y_modal1 = V_eq(2) + epsilon*U_mode1(2)*cos(omega_n1*tlist_nonlinear);
theta_modal1 = V_eq(3) + epsilon*U_mode1(3)*cos(omega_n1*tlist_nonlinear);
plot(x_modal1(1,:), y_modal1(1,:), 'r', 'LineWidth', 3)

figure()
subplot(3,1,1)
plot(tlist_nonlinear, Vlist_nonlinear(1,:))
title('Time vs. X')
subplot(3,1,2)
plot(tlist_nonlinear, Vlist_nonlinear(2,:))
title('Time vs. Y')
subplot(3,1,3)
plot(tlist_nonlinear, Vlist_nonlinear(3,:))
title('Time vs. Theta')

% mode 2
U_mode2 = U_mode(:,2);
omega_n2 = omega_n(2,2);

V0 = V_eq + epsilon*[U_mode2;0;0;0];

%run the integration of nonlinear system
[tlist_nonlinear,Vlist_nonlinear,~,~,~] =...
explicit_RK_variable_step_integration(my_rate_func,tspan,V0,h_ref,DormandPrince,p,error_desired);

figure()
plot(Vlist_nonlinear(1,:), Vlist_nonlinear(2,:), '-b')
hold on
x_modal2 = V_eq(1) + epsilon*U_mode2(1)*cos(omega_n2*tlist_nonlinear);
y_modal2 = V_eq(2) + epsilon*U_mode2(2)*cos(omega_n2*tlist_nonlinear);
theta_modal2 = V_eq(3) + epsilon*U_mode2(3)*cos(omega_n2*tlist_nonlinear);
plot(x_modal2(1,:), y_modal2(1,:), 'r', 'LineWidth', 3)

figure()
subplot(3,1,1)
plot(tlist_nonlinear, Vlist_nonlinear(1,:))
title('Time vs. X')
subplot(3,1,2)
plot(tlist_nonlinear, Vlist_nonlinear(2,:))
title('Time vs. Y')
subplot(3,1,3)
plot(tlist_nonlinear, Vlist_nonlinear(3,:))
title('Time vs. Theta')

% mode 3
U_mode3 = U_mode(:,3);
omega_n3 = omega_n(3,3);

V0 = V_eq + epsilon*[U_mode3;0;0;0];

%run the integration of nonlinear system
[tlist_nonlinear,Vlist_nonlinear,~,~,~] =...
explicit_RK_variable_step_integration(my_rate_func,tspan,V0,h_ref,DormandPrince,p,error_desired);

figure()
plot(Vlist_nonlinear(1,:), Vlist_nonlinear(2,:), '-b')
hold on
x_modal3 = V_eq(1) + epsilon*U_mode3(1)*cos(omega_n3*tlist_nonlinear);
y_modal3 = V_eq(2) + epsilon*U_mode3(2)*cos(omega_n3*tlist_nonlinear);
theta_modal3 = V_eq(3) + epsilon*U_mode3(3)*cos(omega_n3*tlist_nonlinear);
plot(x_modal3(1,:), y_modal3(1,:), 'r', 'LineWidth', 3)

figure()
subplot(3,1,1)
plot(tlist_nonlinear, Vlist_nonlinear(1,:))
title('Time vs. X')
subplot(3,1,2)
plot(tlist_nonlinear, Vlist_nonlinear(2,:))
title('Time vs. Y')
subplot(3,1,3)
plot(tlist_nonlinear, Vlist_nonlinear(3,:))
title('Time vs. Theta')
