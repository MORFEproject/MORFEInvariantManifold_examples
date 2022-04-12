function [bialt_M1,bialt_M2,bialt_M3,bialt_M4]= bialtaa(nphase)
%
% Computes indices of bialternate product of A and A
%
% Syntax: bialtaa(n)

for j=1:(nphase-1)
    i=(j+1):nphase;
    index1=((nphase-1)*nphase/2)-((nphase-j)*(nphase-j+1)/2)+(i-j);
    for q=1:(nphase-1)
        k=(q+1):nphase;
        index2=((nphase-1)*nphase/2)-((nphase-q)*(nphase-q+1)/2)+(k-q);

        k1=repmat(k,size(index1,2),1);
        i1=repmat(i',1,size(index2,2));
        bialt_M1(index1,index2)=(k1-1)*nphase+i1;
        bialt_M2(index1,index2)=(q-1)*nphase+j;
        bialt_M3(index1,index2)=(q-1)*nphase+i1;
        bialt_M4(index1,index2)=(k1-1)*nphase+j;
    end%q
end %j
