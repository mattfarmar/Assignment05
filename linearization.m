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

omega_n = sqrt(omega_n);

t_list = linspace(0,10,500);


% MODE 1
U_mode1 = U_mode(:,1);
omega_n1 = omega_n(1,1);

epsilon = 0.005;
V0 = V_eq + epsilon*[U_mode1;0;0;0];


%run the integration of nonlinear system
[tlist_nonlinear,Vlist_nonlinear,~,~,~] =...
explicit_RK_variable_step_integration(my_rate_func,tspan,V0,h_ref,DormandPrince,p,error_desired);

figure()
plot(Vlist_nonlinear(1,:), Vlist_nonlinear(2,:), '-b')
hold on
x_modal1 = V_eq(1) + epsilon*U_mode1(1)*cos(omega_n1*t_list);
y_modal1 = V_eq(2) + epsilon*U_mode1(2)*cos(omega_n1*t_list);
theta_modal1 = V_eq(3) + epsilon*U_mode1(3)*cos(omega_n1*t_list);
plot(x_modal1(1,:), y_modal1(1,:), 'r', 'LineWidth', 3)

figure()
subplot(3,1,1)
plot(tlist_nonlinear, Vlist_nonlinear(1,:), 'b')
hold on
plot(t_list, x_modal1, 'r')
title('Time vs. X (Mode 1)')
subplot(3,1,2)
plot(tlist_nonlinear, Vlist_nonlinear(2,:), 'b')
hold on
plot(t_list, y_modal1, 'r')
title('Time vs. Y (Mode 1)')
subplot(3,1,3)
plot(tlist_nonlinear, Vlist_nonlinear(3,:), 'b')
hold on
plot(t_list, theta_modal1, 'r')
title('Time vs. Theta (Mode 1)')
legend('Nonlinear', 'Modal')

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
x_modal2 = V_eq(1) + epsilon*U_mode2(1)*cos(omega_n2*t_list);
y_modal2 = V_eq(2) + epsilon*U_mode2(2)*cos(omega_n2*t_list);
theta_modal2 = V_eq(3) + epsilon*U_mode2(3)*cos(omega_n2*t_list);
plot(x_modal2(1,:), y_modal2(1,:), 'r', 'LineWidth', 3)

figure()
subplot(3,1,1)
plot(tlist_nonlinear, Vlist_nonlinear(1,:), 'b')
hold on
plot(t_list, x_modal2, 'r')
title('Time vs. X (Mode 2)')
subplot(3,1,2)
plot(tlist_nonlinear, Vlist_nonlinear(2,:), 'b')
hold on
plot(t_list, y_modal2, 'r')
title('Time vs. Y (Mode 2)')
subplot(3,1,3)
plot(tlist_nonlinear, Vlist_nonlinear(3,:), 'b')
hold on
plot(t_list, theta_modal2, 'r')
title('Time vs. Theta (Mode 2)')
legend('Nonlinear', 'Modal')

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
x_modal3 = V_eq(1) + epsilon*U_mode3(1)*cos(omega_n3*t_list);
y_modal3 = V_eq(2) + epsilon*U_mode3(2)*cos(omega_n3*t_list);
theta_modal3 = V_eq(3) + epsilon*U_mode3(3)*cos(omega_n3*t_list);
plot(x_modal3(1,:), y_modal3(1,:), 'r', 'LineWidth', 3)

figure()
subplot(3,1,1)
plot(tlist_nonlinear, Vlist_nonlinear(1,:), 'b')
hold on
plot(t_list, x_modal3, 'r')
title('Time vs. X (Mode 3)')
subplot(3,1,2)
plot(tlist_nonlinear, Vlist_nonlinear(2,:), 'b')
hold on
plot(t_list, y_modal3, 'r')
title('Time vs. Y (Mode 3)')
subplot(3,1,3)
plot(tlist_nonlinear, Vlist_nonlinear(3,:), 'b')
hold on
plot(t_list, theta_modal3, 'r')
title('Time vs. Theta (Mode 3)')
legend('Nonlinear', 'Modal')

%% FUNCTIONS
%updates spring plotting object so that spring is plotted
%with ends located at points P1 and P2
function update_spring_plot(spring_plot_struct,P1,P2)
    dP = P2-P1;
    R = [dP(1),-dP(2)/norm(dP);dP(2),dP(1)/norm(dP)];
    plot_pts = R*spring_plot_struct.zig_zag;
    set(spring_plot_struct.line_plot,...
    'xdata',plot_pts(1,:)+P1(1),...
    'ydata',plot_pts(2,:)+P1(2));
    set(spring_plot_struct.point_plot,...
    'xdata',[P1(1),P2(1)],...
    'ydata',[P1(2),P2(2)]);
end
%create a struct containing plotting info for a single spring
%INPUTS:
%num_zigs: number of zig zags in spring drawing
%w: width of the spring drawing
function spring_plot_struct = initialize_spring_plot(num_zigs,w)
    spring_plot_struct = struct();
    zig_ending = [.25,.75,1; ...
    -1,1,0];
    zig_zag = zeros(2,3+3*num_zigs);
    zig_zag(:,1) = [-.5;0];
    zig_zag(:,end) = [num_zigs+.5;0];
    for n = 0:(num_zigs-1)
        zig_zag(:,(3+3*n):2+3*(n+1)) = zig_ending + [n,n,n;0,0,0];
    end
    zig_zag(1,:)=(zig_zag(1,:)-zig_zag(1,1))/(zig_zag(1,end)-zig_zag(1,1));
    zig_zag(2,:)=zig_zag(2,:)*w;
    spring_plot_struct.zig_zag = zig_zag;
    spring_plot_struct.line_plot = plot(0,0,'k','linewidth',2);
    spring_plot_struct.point_plot = plot(0,0,'ro','markerfacecolor','r','markersize',7);
end

