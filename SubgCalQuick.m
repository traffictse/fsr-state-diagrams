function SubgCalQuick(degree,option_plot,option_cal)
%calculate all connected subgraphs of a given graph
%output each connected subgraph and distribution table of their sizes
%along with a matrix of conjugate points' distribution
%also plot the given graph and another graph provided by conjmat
%option_plot has four values, namely 1, 2, 3 and 4
%1 is to plot figure1
%2 is to plot figure2
%3 is to plot figure1 and figure2
%4 is not to plot any figure
%option_cal has two values, namely 1 and 2
%1 is to calculate successors with a recursion
%2 is to calculate successors in finite fields with a polynomial
tic;
fid=fopen('datas.txt','w+');

switch option_cal
    case 1
        Successor=SuccCalQuick(degree);
    case 2
        Successor=SuccCalFF(degree);
    otherwise
        error('Invalid input of option_cal! Either 1 or 2!');
end

G1=digraph(Successor(1,:)',Successor(2,:)');

%{
f=fix((sum(Successor)-1)/2).*(fix((sum(Successor)-1)/2)+mod(sum(Successor)-1,2))+min(Successor(1,:),Successor(2,:));
[~,ia,~]=unique(f);
edge=Successor(:,ia);
G2=graph(edge(1,:),edge(2,:));
bins=conncomp(G2);
[Ele,Freq]=CalFreq(bins);
%}
weak_bins=conncomp(G1,'Type','weak');
[Ele,Freq]=CalFreq(weak_bins);
%display(Freq);
%Ele represents mark numbers of all subgraphs
%Freq represents sizes of all subgraphs
[x,y]=CalFreq(Freq');
fprintf(fid,['Distribution Table of Sizes of Connected Subgraphs:' '\r\n']);
fprintf(['Distribution Table of Sizes of Connected Subgraphs:' '\n']);
for j=1:length(x)
    fprintf(fid,[num2str(x(j)) ' ' num2str(y(j)) '\r\n']);
    fprintf([num2str(x(j)) ' ' num2str(y(j)) '\n']);
end
fprintf(fid,['\r\n' 'Totally there are ' num2str(length(Ele)) ' connected subgraphs:' '\r\n']);
fprintf(['\n' 'Totally there are ' num2str(length(Ele)) ' connected subgraphs.' '\n']);
fclose(fid);
%output distribution table of sizes of all connected subgraphs

%Subgraph(1:length(Freq),1:max(Freq))=0;
%for s=1:length(Freq)
%    Subgraph(s,1:Freq(s))=find(weak_bins==Ele(s));
%    %display(find(weak_bins==Ele(s)));
%end

%calculate length of each circle
for s=1:length(Ele)
    Subgraph=[];
    Subgraph(1:2,1:Freq(s))=Successor(1:2,weak_bins==Ele(s));
    %calculate corresponding circle in each connected subgraph respectively
    curr_node=Subgraph(1,1);
    node_vec(1:Freq(s))=0;
    vec_ind=0;
    while(~ismember(curr_node,node_vec))
        vec_ind=vec_ind+1;
        node_vec(vec_ind)=curr_node;
        curr_node=Subgraph(2,Subgraph(1,:)==curr_node);
    end
    [~,rep_ind]=ismember(curr_node,node_vec);
    cir_len=vec_ind+1-rep_ind;
    fid=fopen('datas.txt','a+');
    fprintf(fid,['\r\n' 'Circle in Connected Subgraph ' num2str(s) ' :' ' (' num2str(cir_len) ')' '\r\n']);
    fclose(fid);
    dlmwrite('datas.txt',node_vec(rep_ind:vec_ind),'-append','delimiter','\t','newline','pc');
end

Conjmat(1:length(Freq),1:length(Freq))=0;

%in order to achieve an efficient algorithm, we rewrite the following part
%{
for t=1:length(Freq)
    for u=1:Freq(t)
        if (Subgraph(t,u)>2^(degree-1))
            [m,~]=find(Subgraph==(Subgraph(t,u)-2^(degree-1)));
        else
            [m,~]=find(Subgraph==(Subgraph(t,u)+2^(degree-1)));
        end
        if (m>t)
            Conjmat(t,m)=Conjmat(t,m)+1;
            Conjmat(m,t)=Conjmat(m,t)+1;
        else if (m==t)
            Conjmat(m,m)=Conjmat(m,m)+1/2;
            end
        end
    end   
end
%}

for t=1:2^(degree-1)
    %it's time to make good use of the array weak_bins
    Conjmat(weak_bins(t),weak_bins(t+2^(degree-1)))=Conjmat(weak_bins(t),weak_bins(t+2^(degree-1)))+1;
    Conjmat(weak_bins(t+2^(degree-1)),weak_bins(t))=Conjmat(weak_bins(t+2^(degree-1)),weak_bins(t))+1;
end
Conjmat(1:size(Conjmat,1)+1:end)=Conjmat(1:size(Conjmat,1)+1:end)/2;
[Ele1,Freq1]=CalFreq(diag(Conjmat));

Assocmat(1:length(Freq),1:length(Freq))=0;

for t=1:2:2^degree
    Assocmat(weak_bins(t),weak_bins(t+1))=Assocmat(weak_bins(t),weak_bins(t+1))+1;
    Assocmat(weak_bins(t+1),weak_bins(t))=Assocmat(weak_bins(t+1),weak_bins(t))+1;
end
Assocmat(1:size(Assocmat,1)+1:end)=Assocmat(1:size(Assocmat,1)+1:end)/2;
[Ele2,Freq2]=CalFreq(diag(Assocmat));

%Calculation of triple junctions, selfloops and leaves
[Node,indegree]=CalFreq(Successor(2,:));
%display([Node' indegree']);
TriJunc=Node(indegree==2);
if (numel(TriJunc)==0)
    fid=fopen('datas.txt','a+');
    fprintf(fid,['\r\n' 'There is no triple junction!' '\r\n']);
    fclose(fid);
    disp('There is no triple junction!');
else
    fid=fopen('datas.txt','a+');
    fprintf(fid,['\r\n' 'Triple Junctions:' ' (' num2str(length(TriJunc)) ')' '\r\n']);
    fclose(fid);
    dlmwrite('datas.txt',TriJunc,'-append','delimiter','\t','newline','pc');
    [TJx,TJy]=CalFreq(weak_bins(TriJunc));
    for i=1:length(TJx)
        fid=fopen('datas.txt','a+');
        fprintf(fid,['\r\n' 'Triple Junctions in Connected Subgraph ' num2str(i) ' :' ' (' num2str(TJy(i)) ')' '\r\n']);
        fclose(fid);
        dlmwrite('datas.txt',TriJunc(weak_bins(TriJunc)==TJx(i)),'-append','delimiter','\t','newline','pc');
    end
end
res=Successor(1,:)-Successor(2,:);
SelfLoop=find(res==0);
if (numel(SelfLoop)==0)
    fid=fopen('datas.txt','a+');
    fprintf(fid,['\r\n' 'There is no self loop!' '\r\n']);
    fclose(fid);
    disp('There is no self loop!');
else
    fid=fopen('datas.txt','a+');
    fprintf(fid,['\r\n' 'Self Loops:' ' (' num2str(length(SelfLoop)) ')' '\r\n']);
    fclose(fid);
    dlmwrite('datas.txt',SelfLoop,'-append','delimiter','\t','newline','pc');
    [SLx,SLy]=CalFreq(weak_bins(SelfLoop));
    for i=1:length(SLx)
        fid=fopen('datas.txt','a+');
        fprintf(fid,['\r\n' 'Self loops in Connected Subgraph ' num2str(i) ' :' ' (' num2str(SLy(i)) ')' '\r\n']);
        fclose(fid);
        dlmwrite('datas.txt',SelfLoop(weak_bins(SelfLoop)==SLx(i)),'-append','delimiter','\t','newline','pc');
    end
end
SuccUniq=unique(Successor(2,:));
Leaves=(1:2^degree);
Leaves(SuccUniq)=[];
if (numel(Leaves)==0)
    fid=fopen('datas.txt','a+');
    fprintf(fid,['\r\n' 'There is no leaf!' '\r\n']);
    fclose(fid);
    disp('There is no leaf!');
else
    fid=fopen('datas.txt','a+');
    fprintf(fid,['\r\n' 'Leaves:' ' (' num2str(length(Leaves)) ')' '\r\n']);
    fclose(fid);
    dlmwrite('datas.txt',Leaves,'-append','delimiter','\t','newline','pc');
    [Lx,Ly]=CalFreq(weak_bins(Leaves));
    for i=1:length(Lx)
        fid=fopen('datas.txt','a+');
        fprintf(fid,['\r\n' 'Leaves in Connected Subgraph ' num2str(i) ' :' ' (' num2str(Ly(i)) ')' '\r\n']);
        fclose(fid);
        dlmwrite('datas.txt',Leaves(weak_bins(Leaves)==Lx(i)),'-append','delimiter','\t','newline','pc');
    end
end

disp('Calculation completed.');
toc;

tic;
fid=fopen('datas.txt','a+');
fprintf(fid,['\r\n' 'Conjugated Matrix:' '\r\n']);
fclose(fid);
dlmwrite('datas.txt',Conjmat,'-append','delimiter','\t','newline','pc');
fid=fopen('datas.txt','a+');
fprintf(fid,['\r\n' 'Distribution of Diagonal Elements of Conjugated Matrix:' '\r\n']);
fprintf(['\n' 'Distribution of Diagonal Elements of Conjugated Matrix:' '\n']);
fclose(fid);
dlmwrite('datas.txt',[Ele1' Freq1'],'-append','delimiter','\t','newline','pc');
display([Ele1' Freq1']);
fid=fopen('datas.txt','a+');
fprintf(fid,['\r\n' 'Associated Matrix:' '\r\n']);
fclose(fid);
dlmwrite('datas.txt',Assocmat,'-append','delimiter','\t','newline','pc');
fid=fopen('datas.txt','a+');
fprintf(fid,['\r\n' 'Distribution of Diagonal Elements of Associated Matrix:' '\r\n']);
fprintf(['\n' 'Distribution of Diagonal Elements of Associated Matrix:' '\n']);
fclose(fid);
dlmwrite('datas.txt',[Ele2' Freq2'],'-append','delimiter','\t','newline','pc');
display([Ele2' Freq2']);
disp('Data write-in completed.');
toc;

%plot graphs
switch option_plot
    case 1
        tic;
        disp('Ready to plot...Please wait...');
        nodelabel=cellstr(dec2bin(0:2^degree-1));
        set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w');
        plot(G1,'layout','force','marker','o','markersize',4,'nodelabel',nodelabel);%Or layout -> force/layered
        disp('Ploting completed.');
        toc;
    case 2
        tic;
        disp('Ready to plot...Please wait...');
        figure(1)
        Adj2GraphQuick(Conjmat);
        hold on
        figure(2)
        Adj2GraphQuick(Assocmat);
        disp('Ploting completed.');
        toc;
    case 3
        tic;
        disp('Ready to plot figure1...Please wait...');
        figure(1)
        %subplot(3,1,1);
        nodelabel=cellstr(dec2bin(0:2^degree-1));
        plot(G1,'layout','force','marker','o','markersize',4,'nodelabel',nodelabel);
        hold on
        disp('Ready to plot figure2...Please wait...');
        figure(2)
        %subplot(3,1,2);
        Adj2GraphQuick(Conjmat);
        hold on
        disp('Ready to plot figure3...Please wait...');
        figure(3)
        %subplot(3,1,3);
        Adj2GraphQuick(Assocmat);
        disp('Ploting completed.');
        toc;
    otherwise
        disp('No request for plotting!');
end

tic;
if(Conjmat==Assocmat)
    disp('Conjugated matrix is the same as associated matrix!');
else
    disp('Conjugated matrix and associated matrix are different!');
end
if((Conjmat>2)-zeros(size(Conjmat,1))==zeros(size(Conjmat,1)))
    disp('Elements in conjugated matrix are any of 0, 1 or 2.');
else
    disp('Some elements in conjugated matrix are greater than 2.');
end
if((Assocmat>2)-zeros(size(Assocmat,1))==zeros(size(Assocmat,1)))
    disp('Elements in associated matrix are any of 0, 1 or 2.');
else
    disp('Some elements in associated matrix are greater than 2.');
end
%display(Conjmat);
[row1,~]=find(diag(Conjmat)>0);
%display(row1');
Conjmat(row1,:)=[];
Conjmat(:,row1)=[];
%display(Conjmat);
if(numel(Conjmat)==0)
    disp('Conjugated cofactor is empty.');
elseif((Conjmat>2)-zeros(size(Conjmat,1))==zeros(size(Conjmat,1)))
        disp('Elements in conjugated cofactor are any of 0, 1 or 2.');
    else
        disp('Some elements in conjugated cofactor are greater than 2.');
end
%display(Assocmat);
[row2,~]=find(diag(Assocmat)>0);
%display(row2');
Assocmat(row2,:)=[];
Assocmat(:,row2)=[];
%display(Assocmat);
if(numel(Assocmat)==0)
    disp('Associated cofactor is empty.');
elseif((Assocmat>2)-zeros(size(Assocmat,1))==zeros(size(Assocmat,1)))
        disp('Elements in associated cofactor are any of 0, 1 or 2.');
    else
        disp('Some elements in associated cocfactor are greater than 2.');
end
toc;
disp('Verification completed.');
end