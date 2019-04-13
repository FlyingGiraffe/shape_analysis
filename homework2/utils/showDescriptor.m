function h = showDescriptor(X,T,desc)

h = figure;
patch('vertices',X,'Faces',T,'CData',desc,'FaceColor','interp','edgecolor','none');
%view(22,12)
axis equal;
axis off;
cameratoolbar;
colormap spring;
%colorbar('south','Position',[0.05 0.05 0.9 0.05]);

set(gca,'LooseInset',get(gca,'TightInset')+[0.03 0.03 0.03 0.03]);
set(gcf,'Position',[0,0,512,512]);