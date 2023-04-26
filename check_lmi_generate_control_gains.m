clear all;
global K1 K2 K3 K4 L1 L2 L3 L4
%% Parameters
 lf= 1.035;lr=1.655;
 Vx=30;
 m = 1704.7;Iz=3048.1;mu=0.9;
%%  TS configuration
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
C = [1 0 0 0;0 0 1 0;0 0 0 1];

X = sdpvar(4,4); Y = sdpvar(4,4);
for i =1:1:4
    M{i} = sdpvar(1,4,'full');
    N{i} = sdpvar(4,3,'full');
end
for i =1:1:4
    for j =i:1:4
        if i == j
            P{i,j} = sdpvar(4,4);
            Q{i,j} = sdpvar(4,4);            
        else
            P{i,j} = sdpvar(4,4,'full');
            Q{i,j} = sdpvar(4,4,'full');
        end
    end
end
gamma = sdpvar(1,1);
t1 = diag([-0.5 0.5 -0.5 -0.5]); t2 = 2000;
Con = [X >=0, Y>=1.2*eye(4), gamma>=0]; 
for i = 1:1:4
    Di = gamma*(D{i}*D{i}');
    G1 = A{i}*X + B{i}*M{i} + (A{i}*X + B{i}*M{i})' +  gamma*(D{i}*D{i}') - P{i,i} + 2*t1*X;
    G2 = Y*A{i} + N{i}*C + (Y*A{i} + N{i}*C)' - Q{i,i} + 2*t2*Y;
    Con = [Con, G1 <=0,G2 <=0];
end

for i=1:1:3
    for j=i+1:1:4
        Dij = gamma*(D{i}*D{j}'+D{j}*D{i}');
        G1 = A{i}*X+B{j}*M{i} + A{j}*X + B{i}*M{j} + (A{i}*X+B{j}*M{i} + A{j}*X + B{i}*M{j})' + gamma*(D{i}*D{j}'+D{j}*D{i}') - P{i,j} - P{i,j}'+ 4*t1*X;
        G2 = Y*A{i} +Y*A{j} + N{i}*C + N{j}*C + (Y*A{i} +Y*A{j} + N{i}*C + N{j}*C)' - Q{i,j} - Q{i,j}'+ 4*t2*Y;
        Con = [Con, G1 <=0, G2<=0];
    end
end
F = blkvar;
G = blkvar;
H = [0 0 0 0];
for i =1:1:4
    for j=i:1:4
        F(i,j) = Q{i,j};
        G(i,j) = P{i,j};
    end
    G(i,5) = X*H';
end
G(5,5) = -eye(1,1);
F = sdpvar(F);
G = sdpvar(G);
Con = [Con, F<=0,G<=0];


%Check observability rank(obsv(A,C)) - rank(A)
options_sdp=sdpsettings;
options_sdp.solver='sedumi'; % sedumi sdpt3
% options_sdp.sdpt3.maxit  =50;
options_sdp.sedumi.maxiter  =50;
options_sdp.shift=1e-5; % next two lines: numerical parameters
options_sdp.verbose=1;  % =0: suppress screenoutput, =1: allow screen output

crit = 0;
solpb = solvesdp(Con,[],options_sdp);
[primal,dres] = checkset(Con);


if isequal(sum(primal>0),length(primal)) % see its significance in yalmiperror; checkset(pblmi)
  disp('The problem has been found FEASIBLE');
  
    X_sol = double(X); Y_sol = double(Y);
    K1 = double(M{1})*inv(X_sol);
    K2 = double(M{2})*inv(X_sol); 
    K3 = double(M{3})*inv(X_sol);
    K4 = double(M{4})*inv(X_sol);
    L1 = inv(Y_sol)*double(N{1});
    L2 = inv(Y_sol)*double(N{2}); 
    L3 = inv(Y_sol)*double(N{3});
    L4 = inv(Y_sol)*double(N{4});
    gamma_sol = double(gamma);



else
  disp('The problem has been found IIIIIIIINNNFEASIBLE');  
end