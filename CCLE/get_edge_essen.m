function [ edge_essen ] = get_edge_essen( connectivity,y,alpha,w )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
inds = find(sum(connectivity,2) > 0);
C = connectivity(inds, :);
y = y(inds);

m = size(C,2);
%w = 0.5;


graph_edges= zeros(m,2);
for i=1:m
    nds = find(C(:,i)==1);
    graph_edges(i,1) = nds(1);
    graph_edges(i,2) = nds(2);
end

v_imp = y;

B = C;
M = size(B,2);
N = size(B,1);

G = zeros(M,M);
for u=1:N
    x = find(B(u,:) > 0);
    for i=1:length(x)
        for j=i+1:length(x)
            G(x(i), x(j)) = exp(w*v_imp(u));
            G(x(j), x(i)) =  G(x(i), x(j)) ;
        end
    end
    
end

[~, cc] = graphconncomp(sparse(G), 'Directed', 'false');
large_con = mode(cc);
inds = find(large_con~=cc);
G(:, inds) = [];
G(inds,:) = [];

edge_essen = zeros(M,1);
new_inds = large_con==cc;


GX = diag(1./sum(G))*G;
D = diag(sum(GX,2));
GN = alpha*GX + (1-alpha)*D;


% X = GN;
% for j=1:1
%     X = X*X;
% end
% vec = X(1,:)';


[vv,~] = eigs(GN',1);
vec = vv(:,1);

if(isempty(inds))
    edge_essen(1:end) = vec;
else
    edge_essen(new_inds) = vec;
end
 

end

