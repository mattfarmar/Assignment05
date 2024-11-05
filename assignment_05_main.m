%%% ASSIGNMENT 5 %%%

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
    ];
Plist_world = compute_rbt(x,y,theta,Plist_box);

figure(1);
hold on;
plot(Plist_box(:,1),Plist_box(:,2),'b');
plot(Plist_world(:,1),Plist_world(:,2),'r');
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
