%%%%%%%%-------------------------------------------------------------%%%%%%%%%%
%%%%%%%%          Full COMPLEXITY INVARIANT SET COMPUTATION           %%%%%%%%%%
%%%%%%%%  EXAMPLE: Double Integrator With Parametric Uncertainty     %%%%%%%%%%
%%%%%%%%                   Author: Ankit Gupta                       %%%%%%%%%%
%%%%%%%%           Email: ankit.gupta@chalmers.se                    %%%%%%%%%%
%%%%%%%%               Date: 25-September-2017                       %%%%%%%%%%
%%%%%%%%                  System Description                         %%%%%%%%%%
%%%%%%%% x(k+1)=[1+a11                1+a11]* x(k)+[1]*u(k)+[1]*w(k) %%%%%%%%%%
%%%%%%%%        [  0          1+a22/(1+a11)]       [1]      [1]      %%%%%%%%%%
%%%%%%%%-------------------------------------------------------------%%%%%%%%%%

%%% Install Yalmip toolbox for solving SDP
%%% Install MPT toolbox for plotting the sets

clc; clf; warning off; 
clear all;
addpath('./extras')  
Q=[];

%%%%%%%%%% Volume Maximization/ Minimization %%%%%%%%

optimization = 'max';       %% min/max (min: minimum volume RCI set;  max: maximum volume RCI set)


%%%%%%%%%%%%%% Type of Iteration %%%%%%%%%%%%%%%

no_of_iter = 60;   %% numbers of iterations to perform
tolerance = 10^(-7);   %% Set tolerance

%%%%%%%%%% Bounds on parameters %%%%%%%%%%%%%%%
a11=0.25;a22=0.25;a33=0.25;
    
%%%%%%%%%% System LFT Description %%%%%%%%%%%%%%
%%%%%%  x(k+1)=[A+Bp*Delta*(I-Dp*Delta)^(-1)*Ad]*x(k)+[B+Bp*Delta*(I-Dp*Delta)^(-1)*Bd]*u(k)+[E+Bp*Delta*(I-Dp*Delta)^(-1)*Ed]*w(k);

A  = [1   1;
      0   1];

Ad = [0     1;
      1     1;
      0     1;
      0     1];         

B  = [0;
      1];

Bd = zeros(4,1);

E  = [1;
      0];
  
Ed = zeros(4,1);

Bp = [ 0   1    0   0;
      -1   0    1   1];
  
Dp = [-1   0   1   1;
       0   0   0   0;
       0   0   0   0;
       0   0   0   0];
  
%%%%%%%%%%%%%%   State Constraints (Hx<=d)   %%%%%%%%%%%%%%%%%%

H = 0.2*[eye(2);-eye(2)];
d = [1;1;1;1];
X = Polyhedron(H,d);

%%%%%%%%%%%%%%   Input Constraints (Gu<=h)  %%%%%%%%%%%%%%%%%%%

G = 1/3*[1;-1];
h = [1;1];

%%%%%%%%%%%%%%   Input Disturbance (|Dw|<=v) %%%%%%%%%%%%%%%%%%%

v = 1;
D=[1;-1]/a33;
Wc = a33*[1;-1];  nw = size(Wc,1);

%%%%%%%%%%%%%% Initial Polytope %%%%%%%%%%%%%
% P=eye(2);

% P=[  20  20;
%     -20  0;
%       0 -25];
  
% P=[-0.18   -0.55;
%     0.18   -0.55;
%     0.55   -0.18;
%     0.55    0.18];

V = RandPointsSphere(10,2,2);
P = Polyhedron(V);
P = inv(diag(P.b'))*P.A;
% P= [1.3396   -1.4851;
%     1.6790   -1.0868;
%    -1.8783    0.6870;
%    -1.9503    0.4432;
%    -1.4952   -1.3283;
%     0.8995   -1.7863;
%    -1.3396    1.4851;
%    -1.6790    1.0868;
%     1.8783   -0.6870;
%     1.9503   -0.4432;
%     1.4952    1.3283;
%    -0.8995    1.7863];
% P=rand(10,2);
Theta = Polyhedron([P;-P],ones(2*size(P,1),1));
THT     =   EliminateSymmetricVertices(Theta.V);

%%%%%%%%%%%%%% Parametric Uncertainty Description %%%%%%%%%

syms del1 del2

Delta.structure=diag([del1 del1  del1 del2]);

Delta.bounds=[-a11,a11;
              -a22,a22];

%%%%%%%%%%%%%% Delta matrix computation %%%%%%%%%%%%
vector={};V1=[];
for i=1:size(Delta.bounds,1)
    vector{1,i}=Delta.bounds(i,:);
end
dels = combvec(vector{:});
for i = 1:size(dels,2)
    
    Deltas{i} = double(subs(Delta.structure, {del1, del2}, {dels(:,i)'}));

end

%%%%%%%%%%%%%% RCI set Computation initialization %%%%%%%%%%
no_iter=30;
Wk=eye(size(A,1));Zk=eye(size(A,1));

for i= 1:size(THT,1)*size(P,1)
Xk{i}=eye(size(A,1));
end


for iter=1:no_iter

[W,Xk,N,Zk] = FC_RPI_set(A,Ad,B,Bd,E,Ed,Bp,Dp,Deltas,H,G,Wc,P,THT,iter,1,1,Wk,Xk,Zk,0,1);
Wk=W;
K=N*inv(W);

figure(1)
C = Polyhedron([P*inv(W);-P*inv(W)],ones(2*size(P,1),1));
V1 = [V1;volume(C)];
hold on; 
plot(X,'color','r'); hold on;plot(C,'color','g'); hold on; 
%plot(x'*inv(Zk)*x <= 1);
drawnow
hold on;
X0=W*THT';
X01=X0(:,1);
X02=X0(:,2);
X03=X0(:,3);
X04=X0(:,4);
Xop1=[];Xop2=[];Xop3=[];Xop4=[];
for k=1:20
    Xop1=[Xop1,X01];
    Xop2=[Xop2,X02];
    Xop3=[Xop3,X03];
    Xop4=[Xop4,X04];
    Deltak=Deltas{randi(length(Deltas))};
    wk = a33*(-1 + (1+1)*rand(1,1));
    Y1=(A+Bp*Deltak*inv(eye(4)-Dp*Deltak)*Ad+B*K)*X01+E*wk;
    Y2=(A+Bp*Deltak*inv(eye(4)-Dp*Deltak)*Ad+B*K)*X02+E*wk;
    Y3=(A+Bp*Deltak*inv(eye(4)-Dp*Deltak)*Ad+B*K)*X03+E*wk;
    Y4=(A+Bp*Deltak*inv(eye(4)-Dp*Deltak)*Ad+B*K)*X04+E*wk;
    X01=Y1;
    X02=Y2;
    X03=Y3;
    X04=Y4;
end
plot(Xop1(1,:),Xop1(2,:),'LineWidth',2,'Color','black');drawnow; 
hold on;
plot(Xop2(1,:),Xop2(2,:),'LineWidth',2,'Color','black');drawnow; 
hold on;
plot(Xop3(1,:),Xop3(2,:),'LineWidth',2,'Color','black');drawnow; 
hold on;
plot(Xop4(1,:),Xop4(2,:),'LineWidth',2,'Color','black');drawnow; 
hold on;
end
