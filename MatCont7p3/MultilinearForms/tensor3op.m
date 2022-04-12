function tens3 = tensor3op(T,q1,q2,q3,nphase)
%----------------------------------------------------
% This file computes  T*q1*q2*q3, where T is 3rd derivative 
% of the original vectorfield F at x. 
%----------------------------------------------------
  for b=1:nphase
    S1=T(:,:,:,b);
    BB=tensor2op(S1,q1,q2,nphase);
    S(:,b)=BB;
  end
  tens3=S*q3;  
