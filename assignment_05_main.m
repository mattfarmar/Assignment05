%%% ASSIGNMENT 5 %%%
%% Struct
DormandPrince = make_DP_struct();
%% Test compute_spring_force

k = 10;
l0 = 3;
PA = [0 0];
PB = [3 3];

F = compute_spring_force(k,l0,PA,PB);

%% Test Plist_world
theta = 45;
x = 0;
y = 0;
Plist_box = [
    -5 -5;
    -5 5;
    5 5;
    5 -5;
    -5 -5;
    ]';
Plist_world = compute_rbt(x,y,theta,Plist_box);

figure(1);
hold on;
plot(Plist_box(1,:),Plist_box(2,:),'b');
plot(Plist_world(1,:),Plist_world(2,:),'r');
xlim([-10 10]);
ylim([-10 10]);
axis square;

%% Test compute_accel

box_params = struct();
box_params.m = 5;
box_params.I = 10;
box_params.g = -9.81;
box_params.k_list = [10 10 10 10];
box_params.l0_list = [3 3 3 3];
box_params.P_world = [
0 0 0 0;
0 0 0 0;
];
box_params.P_box = [
3 3 3 3;
3 3 3 3;
];

[ax,ay,atheta] = compute_accel(0,0,0,box_params);

%% testing with mode
% initialize figure

box_params = get_box_params();

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

t_in = 0;
V_in = [1; 1; 0; 0; 0; 0];

my_rate_func = @(V_in) box_rate_func(t_in,V_in,box_params);

V_eq = multi_newton_solver(my_rate_func, V_in, true);

J_approx = approximate_jacobian(my_rate_func, V_eq);
my_rate_func = @(t_in,V_in) box_rate_func(t_in,V_in,box_params);

Q = -J_approx(4:6,1:3);
[U_mode, omega_n] = eig(Q);

omega_n = sqrt(omega_n); 


% MODE
U_mode1 = U_mode(:,1);
omega_n1 = omega_n(1,1);
epsilon = 0.5;
V0 = V_eq + epsilon*[U_mode1;0;0;0];


%     run the integration
[t_list,X_list,~, ~, ~] = explicit_RK_variable_step_integration ...
(my_rate_func,tspan,V0,h_ref,DormandPrince,p,error_desired);

num_zigs = 50;
w = .1;
hold on;
%spring_plot_struct = initialize_spring_plot(num_zigs,w);
axis equal; axis square;
axis(3*[-5,5,-5,5]);

box_plot = plot(0,0,'k','linewidth',2);
spring_1 = initialize_spring_plot(num_zigs,w);
spring_2 = initialize_spring_plot(num_zigs,w);
spring_3 = initialize_spring_plot(num_zigs,w);
spring_4 = initialize_spring_plot(num_zigs,w);

spring_list = {spring_1, spring_2, spring_3, spring_4};

tic;

mypath1 = 'C:\Users\ldao\Downloads\';
fname='mode_1_vibration.avi';
input_fname = [mypath1,fname];

% create a videowriter, which will write frames to the animation file
writerObj = VideoWriter(input_fname);

% must call open before writing any frames
open(writerObj);

fig1 = figure(1);
title("Mode 1 Vibration")
hold on

for i = 1:length(t_list)
    for j = 1:length(box_params.P_world)
        %P2 = [X_list(1,i) - x_dist;X_list(2,i) - y_dist]
        x0 = X_list(1,i);
        y0 = X_list(2,i);
        theta0 = X_list(3,i);
     
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
    current_frame = getframe(fig1);
    writeVideo(writerObj,current_frame)
    pause(0.01);
    
end
close(writerObj);
%% Run Simulation
clf
simulate_box()

function simulate_box()
    box_params = get_box_params();

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
%     run the integration
    [t_list,X_list,~, ~, ~] = explicit_RK_variable_step_integration ...
(my_rate_func,tspan,V0,h_ref,DormandPrince,p,error_desired);


%     figure();
%     subplot(2,1,1)
%     hold on;
%     plot(t_list,X_list(1,:),'r');
%     plot(t_list,X_list(2,:),'b');
% 
%     subplot(2,1,2)
%     plot(t_list,X_list(3,:),'r');

    num_zigs = 50;
    w = .1;
    hold on;
    %spring_plot_struct = initialize_spring_plot(num_zigs,w);
    axis equal; axis square;
    axis(3*[-5,5,-5,5]);
    
    box_plot = plot(0,0,'k','linewidth',2);
    spring_1 = initialize_spring_plot(num_zigs,w);
    spring_2 = initialize_spring_plot(num_zigs,w);
    spring_3 = initialize_spring_plot(num_zigs,w);
    spring_4 = initialize_spring_plot(num_zigs,w);
    
    spring_list = {spring_1, spring_2, spring_3, spring_4};
    
    tic;
    
    for i = 1:length(t_list)
        for j = 1:length(box_params.P_world)
            %P2 = [X_list(1,i) - x_dist;X_list(2,i) - y_dist]
            x0 = X_list(1,i);
            y0 = X_list(2,i);
            theta0 = X_list(3,i);
         
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
end

%% spring plot functions
function spring_plotting_example()
    num_zigs = 5;
    w = .1;
    hold on;
    spring_plot_struct = initialize_spring_plot(num_zigs,w);
    axis equal; axis square;
    axis([-3,3,-3,3]);
    for theta=linspace(0,6*pi,1000)
        P1 = [.5;.5];
        P2 = 2*[cos(theta);sin(theta)];
        update_spring_plot(spring_plot_struct,P1,P2)
        drawnow;
        hold on
    end
end
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

%% RK Variable Step Integration
%Runs numerical integration arbitrary RK method using variable time steps
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)
%BT_struct: a struct that contains the Butcher tableau
% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
%p: how error scales with step size (error = k*h^p)
%error_desired: the desired local truncation error at each step
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration
function [t_list,X_list,h_avg, num_evals, percent_failed] = explicit_RK_variable_step_integration ...
(rate_func_in,tspan,X0,h_ref,BT_struct,p,error_desired)
    num_evals = 0;
    t = tspan(1);
    tf = tspan(2);
    X_list = X0;
    t_list = t; % change num_steps
    XA = X0;

    h = h_ref;

    num_failed_steps = 0;

    num_attempted_steps = 0;
    while t<tf
        num_attempted_steps = num_attempted_steps+1;
        t_next = t+h;

        if t_next>tf
            h= tf-t;
            t_next = tf;
        end

        [XB, num_evals_temp, h_next, redo] = explicit_RK_variable_step...
                (rate_func_in,t,XA,h,BT_struct,p,error_desired);

        num_evals = num_evals+num_evals_temp;
        h = h_next;
        
        if ~redo
            XA = XB;
            t = t_next;
            X_list(:,end+1) = XA;
            t_list(end+1) = t;
        else
            num_failed_steps = num_failed_steps+1;
        end

    end

    h_avg = (tspan(2)-tspan(1))/(length(t_list)-1);
    percent_failed = num_failed_steps/num_attempted_steps;

end

%% RK Variable Step
function [XB, num_evals, h_next, redo] = explicit_RK_variable_step...
(rate_func_in,t,XA,h,BT_struct,p,error_desired)
    alpha = 4; % btwn 1.5 and 10, inclusive
    [XB1, XB2, num_evals] = RK_step_embedded(rate_func_in,t,XA,h,BT_struct); %run 1 step of the solver (on original ts)
    h_next = h*min(0.9*(error_desired/norm(XB1-XB2))^(1/p),alpha); % calculate h_next
    XB = XB1;
    estimated_error = norm(XB1 - XB2); % calculate error
    redo = error_desired<estimated_error;
end
%% RK_step_embedded
%This function computes the value of X at the next time step
%for any arbitrary embedded RK method
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%t: the value of time at the current step
%XA: the value of X(t)
%h: the time increment for a single step i.e. delta_t = t_{n+1} - t_{n}
%BT_struct: a struct that contains the Butcher tableau
% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
%OUTPUTS:
%XB1: the approximate value for X(t+h) using the first row of the Tableau
%XB2: the approximate value for X(t+h) using the second row of the Tableau
%num_evals: A count of the number of times that you called
% rate_func_in when computing the next step
function [XB1, XB2, num_evals] = RK_step_embedded(rate_func_in,t,XA,h,BT_struct)
    k = zeros(length(XA),length(BT_struct.B));
    for i = 1:length(BT_struct.B)
        k(:,i) = rate_func_in(t+BT_struct.C(i)*h, XA+h*(k*BT_struct.A(i,:)'));
    end
    XB1 = XA + h*(k*BT_struct.B(1,:)');
    XB2 = XA + h*(k*BT_struct.B(2,:)');
    num_evals = length(BT_struct.B);
end




