function [V] = RandPointsSphere(m,n,r)

%%% This Function returns '2*m' uniformly distributed points on the 'n' dimensional sphere of radius 'r'
ep = 0.1;Z=[];V=[];

Z = -1 + 2*rand(m,n);

for i=1:m
 
    scale = 1/sqrt(sumsqr(Z(i,:)));
    V = [V;scale*Z(i,:)];
    
end

V = r*[V;-V];