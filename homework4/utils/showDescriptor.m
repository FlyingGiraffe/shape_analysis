function h = showDescriptor(mesh, f, title, camera, scale, h)
if nargin < 3
    title = '';
end

if nargin < 4
    camera = [];
end

if nargin < 5
    scale = [];
end

if nargin < 6
    h = figure;
end

X = mesh.vertices;
T = mesh.triangles;
parentFig = []; %getParentFigure(h);

isAxes = strcmp('axes',get(h,'type'));

if isAxes
    if size(f,1) == mesh.numVertices
        p = patch('vertices',X,'Faces',T,'FaceColor','interp','CData',double(f),'edgecolor','none','parent',h); 
    else
        p = patch('vertices',X,'Faces',T,'FaceColor','flat','CData',double(f),'edgecolor','none','parent',h); 
    end
    axis(h,'equal');
    colorbar('peer',h);
    cameratoolbar;%(parentFig,'Show');
else % h is a figure
    figure(h);
    if size(f,1) == mesh.numVertices
        p = patch('vertices',X,'Faces',T,'FaceColor','interp','CData',double(f),'edgecolor','none'); 
    else
        p = patch('vertices',X,'Faces',T,'FaceColor','flat','CData',double(f),'edgecolor','none'); 
    end
    axis equal;
    colorbar;
    cameratoolbar; cameratoolbar('SetCoordSys','none');
end

lightangle(180,-60)
p.FaceLighting = 'gouraud';
p.AmbientStrength = 0.6;
p.DiffuseStrength = 0.7;
p.SpecularStrength = 1;
p.BackFaceLighting = 'unlit';

set(gca,'LooseInset',get(gca,'TightInset')+[0.03 0.03 0.03 0.03]);
set(gcf,'Position',[0,0,512,512])

if ~isempty(title)
    set(h,'name',title);
end

if nargin >= 4 && isfield(camera, 'up')
    campos(camera.position);
    camva(camera.angle);
    camup(camera.up);
end

if nargin >= 5 && length(scale) == 2
    caxis(scale);
end

if isAxes
    axis(h,'off');
else
    axis off;
end

