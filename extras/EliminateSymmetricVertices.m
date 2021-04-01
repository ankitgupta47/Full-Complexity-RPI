function THT = EliminateSymmetricVertices(V)

THT = [];
tol=1e-9;
N=size(V,1);
Index = [];

for i=1:N
    
    a = repmat(V(i,:),N,1);
    b = a+V;
    
    for j=1:N
        
        
        if (norm(b(j,:),inf) < tol)
            Index = [Index; [i,j]];
        end
        
    end

end

temp=[];

for i=1:size(Index,1)
    
    temp = [temp norm(Index(i,:),inf)];

end

Index = unique(temp);

for i = 1:length(Index)
   
    THT = [THT; V(Index(i),:)];
    
end

end