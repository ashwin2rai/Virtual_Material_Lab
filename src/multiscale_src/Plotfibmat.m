
function Plotfibmat(cmatval,minmat,maxmat,cfibval,minfib,maxfib,fignum)

figure(fignum);
colorleg = colormap(jet);
colornum = @(num,colength) ceil((colength-10)*num+5);

x=0;y=0;r=1;
t = (0:0.01:1)'*2*pi;
x1 = r*sin(t);
y1 = r*cos(t);

rectangle('position',[-2.5 -2.5 5 5],'facecolor',colorleg(colornum...
    (((cmatval-minmat)/(maxmat-minmat)),length(colorleg)),:));
hold on;
fill(x+x1,y+y1,'k','facecolor',colorleg(colornum(...
    ((cfibval-minfib)/(maxfib-minfib)),length(colorleg)),:));

axis([-2.5 2.5 -2.5 2.5]); axis square;
h=colorbar;
h.TickLabels=[{'min'} {''} {''} {''} {''} {''} {''} {''} {''} {''} ...
    {'max'}];
