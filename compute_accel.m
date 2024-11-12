%Computes the linear and angular acceleration of the box
%given its current position and orientation
%INPUTS:
%x: current x position of the box
%y: current y position of the box
%theta: current orientation of the box
%box_params: a struct containing the parameters that describe the system
%Fields:
%box_params.m: mass of the box
%box_params.I: moment of inertia w/respect to centroid
%box_params.g: acceleration due to gravity
%box_params.k_list: list of spring stiffnesses
%box_params.l0_list: list of spring natural lengths
%box_params.P_world: 2 x n list of static mounting
% points for the spring (in the world frame)
%box_params.P_box: 2 x n list of mounting points
% for the spring (in the box frame)
%
%OUTPUTS
%ax: x acceleration of the box
%ay: y acceleration of the box
%atheta: angular acceleration of the box
function [ax,ay,atheta] = compute_accel(x,y,theta,box_params)
    F_total = [0; -box_params.m*box_params.g];
    T_total = 0;
    PC = [x; y];
    PB_list = compute_rbt(x,y,theta,box_params.P_box);
    for i=1:length(box_params.k_list)
        k_temp = box_params.k_list(i);
        l0_temp = box_params.l0_list(i);
        PA_temp = box_params.P_world(:,i);
        PB_temp = PB_list(:,i);
        F_temp = compute_spring_force(k_temp,l0_temp,PA_temp,PB_temp);
        F_total = F_total + F_temp;
        moment_arm = (PB_temp-PC);
        T_temp = moment_arm(1)*F_temp(2)-moment_arm(2)*F_temp(1);
%         T_temp = (PB_temp-PC)'*F_temp;
        T_total = T_total + T_temp;
    end
    linear_accel = F_total/box_params.m;
    ax = linear_accel(1);
    ay = linear_accel(2);
    atheta = T_total / box_params.I;
end
