function [Ele,Freq]= CalFreq(Arr)
%Ele is an array including all unique elements in Arr
%Freq is an array representing the corresponding frequency
%This function is both for a vector or a matrix as input
%Ouput are two arrays
Arr=reshape(Arr,[1,numel(Arr)]);
Arr=sort(Arr);
Arr_diff=[diff(Arr) 1];
Arr_diff=find(Arr_diff);
Ele=Arr(Arr_diff);
Freq=diff([0 Arr_diff]);
end
