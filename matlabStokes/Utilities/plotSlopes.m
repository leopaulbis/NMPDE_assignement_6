function plotSlopes(x,y)

for i=1:length(x)-1
    dx=(x(i+1)-x(i)); dy=(y(i+1)-y(i));
    slope=dy/dx;
    xm=(x(i+1)+x(i))/2; ym=(y(i+1)+y(i))/2;
    text(xm+dx*0,ym+dy*0.1,sprintf('%0.1f',slope),'FontSize',12)
end