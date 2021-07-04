function y=SuccCalQuick(degree)
%only for calculating such kind of recursive formulas
%if degree=n
%then x(n+1)=x(1)+x(n-1)*x(n)
%it can be adapted easily to a new simplified formula though
matrix=(dec2bin(0:(2^degree-1))=='1');
%matrix(:,degree+1)=mod(sum(matrix,2),2);%PSR
%matrix(:,degree+1)=mod(matrix(:,1)+matrix(:,degree),2);
%matrix(:,degree+1)=mod(matrix(:,1)+matrix(:,degree-1).*matrix(:,degree),2);
%matrix(:,degree+1)=mod(1+matrix(:,2)+matrix(:,3)+matrix(:,4)+matrix(:,2).*matrix(:,4)+matrix(:,3).*matrix(:,4)+matrix(:,2).*matrix(:,3).*matrix(:,4),2);
%matrix(:,degree+1)=mod(matrix(:,2)+matrix(:,3)+matrix(:,4)+matrix(:,6)+matrix(:,5).*matrix(:,7).*matrix(:,8),2);
%matrix(:,degree+1)=mod(matrix(:,4)+matrix(:,6)+matrix(:,7)+matrix(:,9),2);
%matrix(:,degree+1)=mod(matrix(:,degree-1).*matrix(:,degree),2);
%matrix(:,degree+1)=mod(matrix(:,2).*matrix(:,3)+matrix(:,4)+matrix(:,5),2);
matrix(:,degree+1)=mod(matrix(:,degree-1)+matrix(:,degree),2);
%matrix(:,degree+1)=mod(1+sum(matrix,2),2);%CSR
%matrix(:,degree+1)=mod(sum(matrix,2)+matrix(:,1),2);
%matrix(:,degree+1)=mod(matrix(:,degree),2);
%matrix(:,degree+1)=mod(matrix(:,1)+matrix(:,1).*matrix(:,2),2);
%matrix(:,degree+1)=mod(matrix(:,2),2);
%matrix(:,degree+1)=mod(matrix(:,1)+matrix(:,degree),2);
%matrix(:,degree+1)=mod(matrix(:,1)+matrix(:,degree)+(matrix(:,2)+1).*(matrix(:,3)+1).*(matrix(:,4)+1),2);
%matrix(:,degree+1)=mod(matrix(:,1)+1,2);
%matrix(:,degree+1)=mod(matrix(:,1)+sum(matrix,2),2);
%matrix(:,degree+1)=mod(matrix(:,degree-1).*matrix(:,degree),2);
%matrix(:,degree+1)=mod(1+matrix(:,degree-1).*matrix(:,degree),2);
%matrix(:,degree+1)=mod(matrix(:,1)+matrix(:,degree-1).*matrix(:,degree),2);
%matrix(:,degree+1)=mod(matrix(:,1),2);
%matrix(:,degree+1)=mod(matrix(:,1).*matrix(:,degree)+matrix(:,degree),2);
%matrix(:,degree+1)=mod(matrix(:,degree-1).*matrix(:,degree)+matrix(:,degree),2);
%matrix(:,degree+1)=mod(matrix(:,degree-1).*matrix(:,degree)+matrix(:,degree),2);
%matrix(:,degree+1)=mod(prod(matrix,2),2);
%matrix(:,degree+1)=mod(sum(matrix(:,1:2:end),2),2);%Sum odds
%matrix(:,degree+1)=mod(sum(matrix(:,2:2:end),2),2);%Sum evens
%matrix(:,degree+1)=mod(matrix(:,1).*matrix(:,2)+matrix(:,degree-1).*matrix(:,degree),2);
%matrix(:,degree+1)=mod(matrix(:,1).*sum(matrix(:,2:end),2)+1,2);
%matrix(:,degree+1)=mod(matrix(:,),2)
%matrix(:,degree+1)=mod(matrix(:,1).*matrix(:,2)+matrix(:,3)+matrix(:,4),2);
%matrix(:,degree+1)=mod(matrix(:1)+matrix(:,3),2);

matrix(:,1)=[];
Successor(1:2^degree,1)=(1:2^degree);
Successor(1:2^degree,2)=1+bin2dec(num2str(matrix));
%{
Successor(1:2^degree,2)=1;
for i=1:degree
Successor(1:2^degree,2)=Successor(1:2^degree,2)+matrix(:,i).*2^(degree-i);
end
%}
y=Successor';
end
