function h = showMesh(X,T)

h = figure;
patch('vertices',X,'Faces',T,'CData',X(:,2),'FaceColor','interp','edgecolor','none');
axis equal;
axis off;
light('Position',[1 0 1],'Style','infinite');
lighting phong;
cameratoolbar;
