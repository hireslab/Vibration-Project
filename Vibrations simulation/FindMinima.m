
% This function scans the array Y and returns its minima

function Z = FindMinima(X,Y)

if length(X(:,1))==1
    X = transpose(X);
end
if length(Y(:,1))==1
    Y = transpose(Y);
end

Z=[];
N = length(Y(:,1));
q=0;
for j=1:1:N-2
    y1 = Y(j,1);
    y2 = Y(j+1,1);
    y3 = Y(j+2,1);
    if y1>y2 && y3>y2
        q=q+1;
        Z(q,1) = X(j+1,1);
        Z(q,2) = Y(j+1,1);
    end
end


