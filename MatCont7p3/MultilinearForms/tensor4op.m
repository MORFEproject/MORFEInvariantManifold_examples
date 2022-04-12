function tens4 = tensor4op(T,q1,q2,q3,q4,nphase)
%----------------------------------------------------
% This file computes  T*q1*q2*q3*q4, where T is 4th derivative 
% of the original vectorfield F at x. 
%----------------------------------------------------
  for b=1:nphase
    S1=T(:,:,:,:,b);
    BB=tensor3op(S1,q1,q2,q3,nphase);
    S(:,b)=BB;
  end
  tens4=S*q4;
