function Adj2GraphQuick(Adj)
%Input an adjacent matrix
%Output a weighted undirected graph with/without minumum spanning tree
Adj=tril(Adj);
[start,terminal,weight]=find(Adj);
[m,~]=size(Adj);
NodeLabel=cell(1,m);
for i=1:m
    NodeLabel{i}=['Node',num2str(i)];
end
G=graph(start,terminal,weight);
p=plot(G,'layout','circle','marker','o','markersize',4,'nodelabel',NodeLabel(1:m),'EdgeLabel',G.Edges.Weight);
p.NodeColor='r';
%[T,~]=minspantree(G);
%highlight(p,T);
end