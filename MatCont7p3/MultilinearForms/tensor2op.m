function tens2 = tensor2op(T,q1,q2,nphase)
%----------------------------------------------------
% This file computes  T*q1*q2, where T is 2nd derivative 
% of the original vectorfield F at x. 
%----------------------------------------------------
  for b=1:nphase
    S(:,b)=T(:,:,b)*q1;
  end
  tens2=S*q2;
