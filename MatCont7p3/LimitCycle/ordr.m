function ial = ordr(msh1, msh2)
ial = zeros(1,length(msh2));
for j=1:length(msh2)
  m = find(msh1>msh2(j));
  if m
    ial(j) = min(m);
  else
    ial(j) = length(msh1);
  end
end
ial = ial-1;