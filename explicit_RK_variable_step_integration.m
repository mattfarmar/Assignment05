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
%     h_test = [h_ref, h_ref];
%     h_list = [];
    X_list = X0;
    t_list = t; % change num_steps
    XA = X0;
%     redo = 1;

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