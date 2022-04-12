


% boundary condition
function bc = bordBVP_R2_bc1()
global lds
bc = [eye(lds.nphase) eye(lds.nphase) ];