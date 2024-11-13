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
plot(tlist_nonlinear, Vlist_nonlinear(1,:), 'b')
hold on
plot(tlist_nonlinear, x_modal1, 'r')
title('Time vs. X')
subplot(3,1,2)
plot(tlist_nonlinear, Vlist_nonlinear(2,:), 'b')
hold on
plot(tlist_nonlinear, y_modal1, 'r')
title('Time vs. Y')
subplot(3,1,3)
plot(tlist_nonlinear, Vlist_nonlinear(3,:), 'b')
hold on
plot(tlist_nonlinear, theta_modal1, 'r')
title('Time vs. Theta')
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
x_modal2 = V_eq(1) + epsilon*U_mode2(1)*cos(omega_n2*tlist_nonlinear);
y_modal2 = V_eq(2) + epsilon*U_mode2(2)*cos(omega_n2*tlist_nonlinear);
theta_modal2 = V_eq(3) + epsilon*U_mode2(3)*cos(omega_n2*tlist_nonlinear);
plot(x_modal2(1,:), y_modal2(1,:), 'r', 'LineWidth', 3)

figure()
subplot(3,1,1)
plot(tlist_nonlinear, Vlist_nonlinear(1,:), 'b')
hold on
plot(tlist_nonlinear, x_modal2, 'r')
title('Time vs. X')
subplot(3,1,2)
plot(tlist_nonlinear, Vlist_nonlinear(2,:), 'b')
hold on
plot(tlist_nonlinear, y_modal2, 'r')
title('Time vs. Y')
subplot(3,1,3)
plot(tlist_nonlinear, Vlist_nonlinear(3,:), 'b')
hold on
plot(tlist_nonlinear, theta_modal2, 'r')
title('Time vs. Theta')
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
x_modal3 = V_eq(1) + epsilon*U_mode3(1)*cos(omega_n3*tlist_nonlinear);
y_modal3 = V_eq(2) + epsilon*U_mode3(2)*cos(omega_n3*tlist_nonlinear);
theta_modal3 = V_eq(3) + epsilon*U_mode3(3)*cos(omega_n3*tlist_nonlinear);
plot(x_modal3(1,:), y_modal3(1,:), 'r', 'LineWidth', 3)

figure()
subplot(3,1,1)
plot(tlist_nonlinear, Vlist_nonlinear(1,:), 'b')
hold on
plot(tlist_nonlinear, x_modal3, 'r')
title('Time vs. X')
subplot(3,1,2)
plot(tlist_nonlinear, Vlist_nonlinear(2,:), 'b')
hold on
plot(tlist_nonlinear, y_modal3, 'r')
title('Time vs. Y')
subplot(3,1,3)
plot(tlist_nonlinear, Vlist_nonlinear(3,:), 'b')
hold on
plot(tlist_nonlinear, theta_modal3, 'r')
title('Time vs. Theta')
legend('Nonlinear', 'Modal')

%% SIMULATE WITH BOX

box_params = get_box_params();
t_in = 0;
V_in = [1; 1; 0; 0; 0; 0];

my_rate_func = @(V_in) box_rate_func(t_in,V_in,box_params);

V_eq = multi_newton_solver(my_rate_func, V_in, true);
J_approx = approximate_jacobian(my_rate_func, V_eq);
Q = -J_approx(4:6,1:3);
[U_mode, omega_n] = eig(Q);

U_mode3 = U_mode(:,3);
omega_n3 = omega_n(3,3);

%load the system parameters into the rate function
%via an anonymous function
my_rate_func = @(t_in,V_in) box_rate_func(t_in,V_in,box_params);
x0 = 0;
y0 = 0;
theta0 = 0;
vx0 = 0;
vy0 = -3;
omega0 = -3;
V0 = [x0;y0;theta0;vx0;vy0;omega0];
tspan = [0, 10];

h_ref = 0.001;
p = 3;
error_desired = 0.0001;

DormandPrince = struct();
DormandPrince.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
DormandPrince.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0;...
5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
DormandPrince.A = [0,0,0,0,0,0,0;
1/5, 0, 0, 0,0,0,0;...
3/40, 9/40, 0, 0, 0, 0,0;...
44/45, -56/15, 32/9, 0, 0, 0,0;...
19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,0;...
9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0;...
35/384, 0, 500/1113, 125/192, -2187/6784, 11/84,0];


num_zigs = 50;
w = .1;
hold on;
axis equal; axis square;
axis(3*[-5,5,-5,5]);

box_plot = plot(0,0,'k','linewidth',2);
spring_1 = initialize_spring_plot(num_zigs,w);
spring_2 = initialize_spring_plot(num_zigs,w);
spring_3 = initialize_spring_plot(num_zigs,w);
spring_4 = initialize_spring_plot(num_zigs,w);

spring_list = {spring_1, spring_2, spring_3, spring_4};

tic;
figure()

for i = 1:length(tlist_nonlinear)
    for j = 1:length(box_params.P_world)
        x0 = Vlist_nonlinear(1,i);
        y0 = Vlist_nonlinear(2,i);
        theta0 = Vlist_nonlinear(3,i);
     
        box_params.P_box;
        Plist_world = compute_rbt(x0,y0,theta0,box_params.P_box);
        Plist_box = compute_rbt(x0,y0,theta0,box_params.boundary_pts);

        set(box_plot,'xdata',Plist_box(1,:),'ydata',Plist_box(2,:));
        P2 = Plist_world(:,j);
        P1 = box_params.P_world(:,j);
        update_spring_plot(spring_list{j},P1,P2)
        
        [~]=toc;
    end
    drawnow;
    pause(0.01);
end

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

