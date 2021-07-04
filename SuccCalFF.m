function y=SuccCalFF(degree)
matrix=(dec2bin(0:(2^degree-1))=='1');
%z=PolyMulFF(2,[1 1],[1 1 1],[1 1 0 1],'power',[1 1 1]);
z=PolyMulFF(2,[1 1 1],'power',3);
z=z(1:end-1);
display(z);
matrix(:,degree+1)=mod(sum(matrix(:,z>0),2),2);
matrix(:,1)=[];
Successor(1:2^degree,1)=(1:2^degree);
Successor(1:2^degree,2)=1;
for i=1:degree
Successor(1:2^degree,2)=Successor(1:2^degree,2)+matrix(:,i).*2^(degree-i);
end
y=Successor';
end

