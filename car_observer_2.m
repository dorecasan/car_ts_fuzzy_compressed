function [sys, x0, str, ts] = car_observer_2(t,x,u,flag)
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
    sizes.NumContStates = 4;
    sizes.NumDiscStates = 0;
    sizes.NumOutputs = 4;
    sizes.NumInputs = 5 ;
    sizes.DirFeedthrough = 0;
    sizes.NumSampleTimes = 0;
    sys = simsizes(sizes);
    x0 = [0; 0 ;0; 0];
    str = [];
    ts = [];
end

function sys = mdlDerivatives(t,x,u) 
 lf= 1.035;lr=1.655;
 global L1 L2 L3 L4
 %% TS Configuration
 m = 1704.7;Iz=3048.1;Vx=30;mu=0.9;
 Cf_min = 100000; Cf_max = 110000; Cfs = [Cf_min, Cf_max];
 Cr_min = 70000; Cr_max = 80000; Crs = [Cr_min,Cr_max];
for i = 1:1:2
    Cf = Cfs(i);
    for k=1:1:2
        Cr = Crs(k);
        a11 = -2*mu*(Cf+Cr)/(m*Vx);a21=2*mu*(lr*Cr-lf*Cf)/(Iz*Vx);
        a12 = 2*mu*(Cf+Cr)/m; a22 = 2*mu*(lf*Cf-lr*Cr)/Iz;
        a13 = 2*mu*(lr*Cr-lf*Cf)/(m*Vx);a23=-2*mu*(Cf*lf^2+Cr*lr^2)/(Iz*Vx);
        b1= 2*mu*Cf/m;b2=2*mu*lf*Cf/Iz;
        d1=2*mu*(lr*Cr-lf*Cf)/(m*Vx)-Vx;d2=-2*mu*(Cf*lf^2+Cr*lr^2)/(Iz*Vx);
        Ai = [0 1 0 0;
               0 a11 a12 a13;
               0 0 0 1;
               0 a21 a22 a23];
        Bi = [0; b1;0;b2];
        Di = [0;d1;0;d2];
        A{k+2*(i-1)} = Ai;
        B{k+2*(i-1)} = Bi;
        D{k+2*(i-1)} = Di;  
    end
end
C = [1 0 0 0;0 0 1 0; 0 0 0 1];
Cf=105000;Cr=75000;
%% Membership function
fm_max = (Cf-Cf_min)/(Cf_max-Cf_min);fm_min = (Cf_max-Cf)/(Cf_max-Cf_min);
fz_max = (Cr-Cr_min)/(Cr_max-Cr_min);fz_min = (Cr_max-Cr)/(Cr_max-Cr_min);
w1 = fm_min*fz_min;w2 = fm_min*fz_max;
w3 = fm_max*fz_min;w4 = fm_max*fz_max;
%% Observer gain
L{1} = L1;L{2} = L2;L{3}=L3;L{4}=L4;
x_estimated = [x(1);x(2);x(3);x(4)];
control = u(1);
y = [u(2);u(3);u(4)];
dw_des = u(5);
y_hat = C*x_estimated;
%-------------------------------------------------------------------------
for i=1:1:4
 Q{i} = A{i}*x_estimated + B{i}*control - L{i}*(y-y_hat) + D{i}*dw_des;
end
dx_estimated = (w1*Q{1}+w2*Q{2}+w3*Q{3}+w4*Q{4})/(w1+w2+w3+w4);

sys(1) = dx_estimated(1);
sys(2) = dx_estimated(2);
sys(3) = dx_estimated(3);
sys(4) = dx_estimated(4);

end

function sys = mdlOutputs(t,x,u)
sys(1)=x(1);
sys(2)=x(2);
sys(3)=x(3);
sys(4)=x(4);  
end

