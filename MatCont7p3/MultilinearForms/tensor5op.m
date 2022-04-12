function tens5 = tensor5op(T,q1,q2,q3,q4,q5,nphase)
%----------------------------------------------------
% This file computes  T*q1*q2*q3*q4*q5, where T is 5th derivative 
% of the original vectorfield F at x. 
%----------------------------------------------------
  for b=1:nphase
    S1=T(:,:,:,:,:,b);
    BB=tensor4op(S1,q1,q2,q3,q4,nphase);
    S(:,b)=BB;
  end
  tens5=S*q5;
  
