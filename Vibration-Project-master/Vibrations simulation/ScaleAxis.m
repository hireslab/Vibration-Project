
function ScaleAxis(s)

axis tight
Ylim = ylim;
Xlim = xlim;
dx = abs(Xlim(2)-Xlim(1));
dy = abs(Ylim(2)-Ylim(1));

set(gca,'Xlim',[Xlim(1)-s*dx,Xlim(2)+s*dx])
set(gca,'Ylim',[Ylim(1)-s*dy,Ylim(2)+s*dy])


