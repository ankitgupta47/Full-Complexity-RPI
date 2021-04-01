function [W,X1,N,Z] = FC_RPI_set(A,Ad,B,Bd,E,Ed,Bp,Dp,Deltas,H,G,Wc,P,THT,iter,sai1,sai2,Wk,Xk,Zk,ellipse_cost,max)

%%%%%%%%-------------------------------------------------------------%%%%%%%%%%
%%%%%%%%               Full Complexity RCI set (u = Kx)              %%%%%%%%%%
%%%%%%%%                   Author: Ankit Gupta                       %%%%%%%%%%
%%%%%%%%           Email: ankit.gupta@chalmers.se                    %%%%%%%%%%
%%%%%%%%                   Date: 21-June-2019                        %%%%%%%%%%
%%%%%%%%-------------------------------------------------------------%%%%%%%%%%    

%%%  x(k+1) = A(Delta)x(k)+B(Delta)u(k)+ E(Delta)w(k),
%%% where , [A(Delta)| B(Delta)| E(Delta)] = [A| B| E]+Bp*Delta*(I-Dp*Delta)[Ad| Bd| Ed]
%%%
%%%  Deltas = cell containing all combination of Delta at the vertex.
%%%  H = state constraint matrix of form Hx<=1
%%%  G = input constraint matrix of form Gu<=1
%%%  Wc = vertices of the input disturbance
%%%  P = initial guess of the invariant set
%%%  THT = vertices of the initial guess of the invariant set
%%%  iter = iteration number in the for-loop
%%%  sai1 = parameter (in invariance condition) for line search to find feasible initial solution (default 1)
%%%  sai2 = parameter (in ellipse (max) condition) for line search to find feasible initial solution (default 1)
%%%  Wk = previous value of 'W' (if iter=1 then set Wk=I)
%%%  Xk = previous value of 'X' (if iter=1 then set Xk's=I)
%%%  Zk = previous value of 'Z' (if iter=1 then set Zk=I)
%%%  ellipse_cost = set to '1' is ellipse is used for volume maximization/minimization else set it to '0'
%%%  max = set to '1' if for volume maximization, for minimization set it '0'

%%%  output
%%% W = current rotation and scaling matrix for the polytope
%%% X1 = for the next iteration
%%% N = Needed to construct control law K = N*inv(W)
%%% Z = slack variable for optimization

ep=10^-8;   %%% tolerance
Cons = [];
%-------------------------------------------------------------------------%
%%% extracting the index %%%
%-------------------------------------------------------------------------%
[nx,nu] =   size(B);
nh      =   size(H,1);
ng      =   size(G,1);
Nu      =   size(THT,1);
p       =   size(P,1);
ei      =   eye(p);
np      =   size(Bp,2);
nw      =   size(Wc,1);

%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
%%% Creating variables %%%
%-------------------------------------------------------------------------%

phi = diag(sdpvar(p,1,'full')); 
Cons = [Cons, phi>=ep*eye(p)];

mu1 = diag(sdpvar(p,1,'full')); 
Cons = [Cons, mu1>=ep*eye(p)];

W = sdpvar(nx,nx,'full');
N = sdpvar(nu,nx,'full');
Z = sdpvar(nx,nx);

for i = 1:p

    X{i} = sdpvar(nx,nx);
    Cons = [Cons, X{i}>=ep*eye(nx)];
    
    
    V{i} = sdpvar(nx,nx,'full');
    
    Di{i} = diag(sdpvar(p,1,'full'));
    Cons = [Cons, Di{i}>=ep*eye(p)];
    
    
    Q{i} = sdpvar(np,np);
    R{i} = sdpvar(np,np);
    S{i} = sdpvar(np,np,'full');
    M{i} = [Q{i} S{i};S{i}' R{i}];
    
    
    if(iter==1)   
        l{i}=sai1*eye(nx)*W+W'*sai1*eye(nx)-sai1*eye(nx)*X{i}*sai1*eye(nx); 
        
        if(ellipse_cost==1  && max==1)
            m1=sai2*eye(nx)*W+W'*sai2*eye(nx)-sai2*eye(nx)*Z*sai2*eye(nx);
        end
    
    else
        
        l{i}=Wk'*inv(Xk{i})'*W+W'*inv(Xk{i})*Wk-Wk'*inv(Xk{i})'*X{i}*inv(Xk{i})*Wk; 
        if(ellipse_cost==1  && max==1)
            m1=Wk'*inv(Zk)'*W+W'*inv(Zk)*Wk-Wk'*inv(Zk)'*Z*inv(Zk)*Wk;
        end
    end
    
end
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
%%% Polya relaxation based condition %%%
%-------------------------------------------------------------------------%

for j = 1:p 
    
    for i = 1:length(Deltas)
        
        Cons = [Cons, [Deltas{i}';eye(size(Bp,2))]'*M{j}*[Deltas{i}';eye(size(Bp,2))]<=0];
        
        for k=i+1:length(Deltas)
            
            Cons = [Cons, [Deltas{i}';eye(size(Bp,2))]'*M{j}*[Deltas{k}';eye(size(Bp,2))]+[Deltas{k}';eye(size(Bp,2))]'*M{j}*[Deltas{i}';eye(size(Bp,2))]<=0];
            
        end
    end
end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
%%% State and Input Constaraints, and Invariance Conditions
%-------------------------------------------------------------------------%

for j = 1:Nu
    
    %%% State and Input Constaraints
  
    Cons = [Cons, [H*W*THT(j,:)'<= ones(nh,1), G*N*THT(j,:)' <= ones(ng,1)]];
    
    
end

for i = 1:p
    
     %%% Condition for invariance %%%
            zz1 = [l{i}                       phi(i,i)*P'*ei(:,i);
                    (phi(i,i)*P'*ei(:,i))'               phi(i,i)];
                
            Cons = [Cons, zz1>=ep*eye(nx+1)];
  
    for k = 1:nw
        
       b = ones(p,1);
       zz2 = [phi(i,i)-b'*Di{i}*b     zeros(1,nx)      (E*Wc(k,:)')'          zeros(1,nx)                                   (Ed*Wc(k,:)')';
                 zeros(1,nx)'         P'*Di{i}*P         (A*W+B*N)'          zeros(nx,nx)                                     (Ad*W+Bd*N)';
                  E*Wc(k,:)'           (A*W+B*N)    V{i}+V{i}'+Bp*R{i}*Bp'       V{i}'                                Bp*S{i}'+Bp*R{i}*Dp';
                 zeros(1,nx)'        zeros(nx,nx)'           V{i}                X{i}                                        zeros(np,nx)';
                  Ed*Wc(k,:)'        (Ad*W+Bd*N)       S{i}*Bp'+Dp*R{i}*Bp'  zeros(np,nx)               Q{i}+Dp*S{i}'+S{i}*Dp'+Dp*R{i}*Dp'];
       
        
                
            Cons = [Cons, zz2>=0];
    end
    
end

%-------------------------------------------------------------------------%

% Sol_opt = [];
Sol_opt = sdpsettings('solver', 'sedumi');


%-------------------------------------------------------------------------%
%%% Cost Function and solving SDP %%%
%-------------------------------------------------------------------------%
if (max==1)

    %%% Volume Maximization %%%
    
    if(ellipse_cost==0)
        
        if(iter==1)
        
            optimize(Cons,-logdet(W+W'),Sol_opt);
        
        else
            
            Cons = [Cons, [W'*Wk+Wk'*W-Wk'*Wk>=Z,Z>=0]];
            optimize(Cons,-logdet(Z),Sol_opt);
        
        end
        
        
    else
        
        for j=1:p
            
            
            xx1=[2*mu1(j,j)-1           (mu1(j,j)*P'*ei(:,j))'; 
                 mu1(j,j)*P'*ei(:,j)                        m1];
         
            Cons = [Cons, xx1>=0];
        
        end
        
        optimize(Cons,-logdet(Z),Sol_opt);
    
    end
    
else
    
    %%% Volume Minimization %%%
    
    if(ellipse_cost==0)
        
        Cons = [Cons, [Z  W';W  eye(nx)]>=0];
        optimize(Cons, trace(Z),Sol_opt);
   
    else
        
        for j=1:Nu
            
            xx1=[1           (W*THT(j,:)')'; 
                 W*THT(j,:)'              Z];
             
            Cons = [Cons, xx1>=0];
        
        end
        
        optimize(Cons,trace(Z),Sol_opt);
    
    end
end
%-------------------------------------------------------------------------%


if(ellipse_cost == 1)
    Z=value(Z);
else 
    Z=eye(nx);
end

W = value(W);
for i=1:length(X) 
    X1{i}=value(X{i});
end
N=value(N);