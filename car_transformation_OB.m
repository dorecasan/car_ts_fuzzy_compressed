function [sys, x0, str, ts] = car_transformation_OB(t,x,u,flag)
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
    sizes.NumContStates = 1;
    sizes.NumDiscStates = 0;
    sizes.NumOutputs = 4;
    sizes.NumInputs = 6 ;
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 0;
    sys = simsizes(sizes);
    x0 = [0];
    str = [];
    ts = [];
end

function sys = mdlDerivatives(t,x,u) 
%% Parameters
%% Inputs
yaw = u(2); w_des = u(5); 
%% Model
e = yaw - w_des;
sys(1) = e;

end

function sys = mdlOutputs(t,x,u)
 y = u(1); yaw = u(2); dyaw = u(3); dw_des = u(4); w_des = u(5); dy = u(6);
 Vx = 30;
 sys(1) = y + Vx*x(1);
 sys(2) = dy +Vx*(yaw- w_des);
 sys(3) = yaw - w_des;
 sys(4) = dyaw - dw_des;
end

