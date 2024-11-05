%%% ASSIGNMENT 5 %%%

%% Run Simulation

function simulate_box()
    %define system parameters
    box_params = struct();
    box_params.m = 5;
    box_params.I = 10;
    box_params.k_list = [];
    box_params.l0_list = [];
    box_params.P_world = [];
    box_params.P_box = [];
    %load the system parameters into the rate function
    %via an anonymous function
    my_rate_func = @(t_in,V_in) box_rate_func(t_in,V_in,box_params);
    x0 = 0;
    y0 = 0;
    theta0 = 0;
    vx0 = 0;
    vy0 = 0;
    omega0 = 0;
    V0 = [x0;y0;theta0;vx0;vy0;omega0];
    tspan = linspace(1,10,100);
    %run the integration
    % [tlist,Vlist] = your_integrator(my_rate_func,tspan,V0,...)
end

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
