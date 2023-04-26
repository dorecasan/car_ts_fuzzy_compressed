function [sys, x0, str, ts] = car_setpoint_OB(t,x,u,flag)
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
    sizes.NumOutputs = 3;
    sizes.NumInputs = 0 ;
    sizes.DirFeedthrough = 0;
    sizes.NumSampleTimes = 1;
    sys = simsizes(sizes);
    x0 = [0];
    str = [];
    ts = [0 0];
    global coef1 coef2 ti tf ti0 Vx tc w2;
    T = 5; ti = 0.5 ; tf =ti+T; tc=2.5; Vx = 30;
    w1 = 0; w2 = 3.75;
    A = [ti^5 ti^4 ti^3 ti^2 ti 1;
     5*ti^4 4*ti^3 3*ti^2 2*ti 1 0;
     20*ti^3 12*ti^2 6*ti 2 0 0;
     tf^5 tf^4 tf^3 tf^2 tf 1;
     5*tf^4 4*tf^3 3*tf^2 2*tf 1 0;
     20*tf^3 12*tf^2 6*tf 2 0 0];

    B = [w1;0;0;w2;0;0];
    coef1  = inv(A)*B;

    B = [w2;0;0;w1;0;0];
    ti0 = ti;
    ti = tf+tc;
    tf = ti+T;
    A = [ti^5 ti^4 ti^3 ti^2 ti 1;
     5*ti^4 4*ti^3 3*ti^2 2*ti 1 0;
     20*ti^3 12*ti^2 6*ti 2 0 0;
     tf^5 tf^4 tf^3 tf^2 tf 1;
     5*tf^4 4*tf^3 3*tf^2 2*tf 1 0;
     20*tf^3 12*tf^2 6*tf 2 0 0];
    coef2  = inv(A)*B;

end
function sys = mdlDerivatives(t,x,u) 
global coef1 coef2 ti tf ti0 Vx tc w2;
if t<ti0
    y = 0; dy =0; ddy = 0;
elseif (ti0<=t) && (t<(ti-tc))
    y = coef1(1)*t^5+coef1(2)*t^4+coef1(3)*t^3+coef1(4)*t^2+coef1(5)*t+coef1(6);
    dy = 5*coef1(1)*t^4+4*coef1(2)*t^3+3*coef1(3)*t^2+2*coef1(4)*t+coef1(5);
    ddy = 20*coef1(1)*t^3+12*coef1(2)*t^2+6*coef1(3)*t+2*coef1(4);
elseif ((ti-tc)<=t) && (t<ti)
    y = w2; dy =0; ddy = 0;
elseif (ti<=t) && (t<tf)
    y = coef2(1)*t^5+coef2(2)*t^4+coef2(3)*t^3+coef2(4)*t^2+coef2(5)*t+coef2(6);
    dy = 5*coef2(1)*t^4+4*coef2(2)*t^3+3*coef2(3)*t^2+2*coef2(4)*t+coef2(5);
    ddy = 20*coef2(1)*t^3+12*coef2(2)*t^2+6*coef2(3)*t+2*coef2(4);
else
    y = 0; dy =0; ddy = 0;
end
psi = Vx*Vx*ddy/((Vx^2+dy^2)^(3/2));
sys(1) = psi;
end
function sys = mdlOutputs(t,x,u)
global coef1 coef2 ti tf ti0 Vx tc w2;
if t<ti0
    y = 0; dy =0; ddy = 0;
elseif (ti0<=t) && (t<(ti-tc))
    y = coef1(1)*t^5+coef1(2)*t^4+coef1(3)*t^3+coef1(4)*t^2+coef1(5)*t+coef1(6);
    dy = 5*coef1(1)*t^4+4*coef1(2)*t^3+3*coef1(3)*t^2+2*coef1(4)*t+coef1(5);
    ddy = 20*coef1(1)*t^3+12*coef1(2)*t^2+6*coef1(3)*t+2*coef1(4);
elseif ((ti-tc)<=t) && (t<ti)
    y = w2; dy =0; ddy = 0;
elseif (ti<=t) && (t<tf)
    y = coef2(1)*t^5+coef2(2)*t^4+coef2(3)*t^3+coef2(4)*t^2+coef2(5)*t+coef2(6);
    dy = 5*coef2(1)*t^4+4*coef2(2)*t^3+3*coef2(3)*t^2+2*coef2(4)*t+coef2(5);
    ddy = 20*coef2(1)*t^3+12*coef2(2)*t^2+6*coef2(3)*t+2*coef2(4);
else
    y = 0; dy =0; ddy = 0;
end
psi = Vx*Vx*ddy/((Vx^2+dy^2)^(3/2));
sys(1) = psi;
sys(2) = x(1);
sys(3) = y;
end

