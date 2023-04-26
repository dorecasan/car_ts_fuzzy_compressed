function [sys, x0, str, ts] = car_ts_controller_OB(t,x,u,flag)
    switch flag
    case 0
        [sys, x0, str,ts] = mdlInitializeSizes;
    case 1
        sys = mdlDerivatives(t,x,u);
    case 3
        sys = mdlOutputs(t,x,u);
    case {2,4,9}
        sys = [];
    otherwise
        error(['Unhandled flag = ', num2str(flag)]);
    end
end
function [sys,x0,str,ts]=mdlInitializeSizes
    sizes = simsizes;
    sizes.NumContStates = 0;
    sizes.NumDiscStates = 0;
    sizes.NumOutputs = 2;
    sizes.NumInputs =5 ;
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 0;
    sys = simsizes(sizes);
    x0 = [];
    str = [];
    ts = [];
    global u_prev
    u_prev = 0;
end


function sys = mdlOutputs(t,x,u) 
%% Parameters
Cf_min = 100000; Cf_max = 110000; 
Cr_min = 70000; Cr_max = 80000; 
global K1 K2 K3 K4 u_prev
Cf=105000;Cr=75000;
%% Membership function
fm_max = (Cf-Cf_min)/(Cf_max-Cf_min);fm_min = (Cf_max-Cf)/(Cf_max-Cf_min);
fz_max = (Cr-Cr_min)/(Cr_max-Cr_min);fz_min = (Cr_max-Cr)/(Cr_max-Cr_min);
w1 = fm_min*fz_min;w2 = fm_min*fz_max;
w3 = fm_max*fz_min;w4 = fm_max*fz_max;
%% Control gains

Kf = (w1*K1+w2*K2+w3*K3+w4*K4)/(w1+w2+w3+w4);
%% Variables
e1 = u(1); de1 = u(2) ; e2 =u(3); de2 = u(4);dw_des = u(5);
x1 = e1; x2 = de1 ; x3 = e2; x4 = de2;
q = [x1;x2;x3;x4];
%% Control signal
u_hat = Kf*q;
control = u_hat;   
% a_max = 0.005;
% if abs(control - u_prev) > a_max
%     control = sign(control-u_prev)*a_max +u_prev;
% end
% u_prev = control;
sys(1)  = control;
sys(2)  = control; 
 
end

