function myarrow_v(x,y)
ax = gca;
axpos = get(ax, 'Position')
X = get(gca,'XLim')
Y = get(gca,'YLim') 
difX = X(2) - X(1);
difY = Y(2) - Y(1);
newx = x./difX;
newy = y./difY;
HLv=7;
HWv=4;
annotation('arrow',[newx(1)*axpos(3)+axpos(1) newx(2)*axpos(3)+axpos(1)],[newy(1)*axpos(4)+axpos(2) newy(2)*axpos(4)+axpos(2)],'HeadLength',HLv,'HeadWidth',HWv);
end