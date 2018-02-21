
% This function generates a row vector from a to b, having length npoints

function v = vec(a,b,npoints)

dn = (b-a)/(npoints-1);
v = a:dn:b;


