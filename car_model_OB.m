function [sys, x0, str, ts] = car_model_OB(t,x,u,flag)
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
    sizes.NumContStates = 5;
    sizes.NumDiscStates = 0;
    sizes.NumOutputs = 5;
    sizes.NumInputs =1 ;
    sizes.DirFeedthrough = 0;
    sizes.NumSampleTimes = 0;
    sys = simsizes(sizes);
    x0 = [0 0 0 0 0];
    str = [];
    ts = [];
end

function sys = mdlDerivatives(t,x,u) 
%% Parameters
 m = 1704.7;Iz=3048.1;lf= 1.035;lr=1.655;
 Cf=105850;Cr=79030;mu=0.9;Vx=30;

 %% Model elements
 a11 = -2*mu*(Cf+Cr)/(m*Vx);a21=2*mu*(lr*Cr-lf*Cf)/(Iz*Vx);
 b1= 2*mu*Cf/m;b2=2*mu*lf*Cf/Iz;
 d1=2*mu*(lr*Cr-lf*Cf)/(m*Vx)-Vx;d2=-2*mu*(Cf*lf^2+Cr*lr^2)/(Iz*Vx);
 A = [0 1 0 0;
     0 a11 0 d1;
     0 0 0 1;
     0 a21 0 d2];
 B = [0; b1;0;b2];
 %% Variables
 theta = u(1);
 y = x(1);
 dy = x(2);
 yaw = x(3);
 dyaw = x(4);
 %% Model
 q = [y dy yaw dyaw]';
 dq = A*q+B*theta;
 sys(1) = dq(1);
 sys(2) = dq(2);
 sys(3) = dq(3);
 sys(4) = dq(4);
 sys(5) = dy+Vx*yaw;
 end

function sys = mdlOutputs(t,x,u)
 sys(1) = x(1);
 sys(2) = x(2);
 sys(3) = x(3);
 sys(4) = x(4);
 sys(5) = x(5);
end

