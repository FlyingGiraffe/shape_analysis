clear; clc; close all;
%% Problem 3.a
% Sample a Lissajous curve
a = 3; b = 2;
delta = pi/2;
n = 100;

t = linspace(0,2*pi,n);
x = sin(a*t*delta);
y = sin(b*t);
z = .5*t;
xyz = [x; y; z];

% Compute point differences and midpoints
diff = xyz(:,2:end)-xyz(:,1:end-1);
midpoint = (xyz(:,1:end-1)+xyz(:,2:end))/2;

% binormal = zeros(3,n-2);
%%% YOUR CODE TO COMPUTE THE BINORMAL HERE %%%
diff_len = sqrt(sum(diff.*diff));
T = diff./diff_len; % tangent = normalize(diff)
%binormal = cross(T(:,1:end-1),T(:,2:end));
%binormal = binormal./sqrt(sum(binormal.*binormal));
binormal = 2*cross(diff(:,1:end-1),diff(:,2:end))./((diff_len(:,1:end-1)+diff_len(:,2:end))+dot(diff(:,1:end-1),diff(:,2:end)));
%%% END HOMEWORK PROBLEM %%%

f1 = figure;
f1.GraphicsSmoothing = 'on';
f1.Renderer = 'painters';
plot3(x,y,z,'linewidth',2,'color','black'); hold on;
quiver3(x(2:end-1),y(2:end-1),z(2:end-1),...
    binormal(1,:),binormal(2,:),binormal(3,:),...
    3,'color','red','linewidth',1);
axis equal;
set(gcf,'Position',[0,0,512,512]);
saveas(gcf,'3a_binormal.png');

%% Problem 3.b
u = zeros(3,n-1);
v = zeros(3,n-1);

h = figure;
h.GraphicsSmoothing = 'on';
h.Renderer = 'painters';
figure(h);
scale=0.25;
plot3(x,y,z,'linewidth',2,'color','black'); hold on;
quiv1 = quiver3(midpoint(1,:),midpoint(2,:),midpoint(3,:),...
    scale*u(1,:),scale*u(2,:),scale*u(3,:),...
    'color','red','linewidth',1);
quiv2 = quiver3(midpoint(1,:),midpoint(2,:),midpoint(3,:),...
    scale*v(1,:),scale*v(2,:),scale*v(3,:),...
    'color','blue','linewidth',1);
hold off;
axis manual;
axis equal;
axis([-2 2 -2 2 0 4]);
set (gcf,'Position',[0,0,512,512]);

writerObj = VideoWriter('3b_bishop_frame.avi');
open(writerObj);
frame = getframe;
for a=0:5:360
    t = diff(:,1)/norm(diff(:,1));
    frame0 = cross(t,[0;0;1]);
    frame1 = cross(frame0,t);
    
    theta = a*pi/180;
    u(:,1) = cos(theta).*frame0+sin(theta).*frame1;
    v(:,1) = -sin(theta).*frame0+cos(theta).*frame1;
    for j=2:(n-1)
        %%% YOUR CODE TO UPDATE u, v HERE %%%
        R = vrrotvec2mat(vrrotvec(T(:,j-1),T(:,j)));
        u(:,j) = R*u(:,j-1);
        v(:,j) = cross(u(:,j),T(:,j));
        %%% END HOMEWORK PROBLEM %%%
    end
    quiv1.UData = scale*u(1,:); quiv1.VData = scale*u(2,:); quiv1.WData = scale*u(3,:); 
    quiv2.UData = scale*v(1,:); quiv2.VData = scale*v(2,:); quiv2.WData = scale*v(3,:);
    drawnow;
    %if mod(a,90) == 0 && ~(a == 360) 
    %    saveas(gcf,['3b_bishop_theta-' int2str(a) '.png']);
    %end
    frame = getframe;
    writeVideo(writerObj,frame);
end
close(writerObj);

%% Problem 3.c
thetas = linspace(0,3*pi,n-1);
m1 = cos(thetas).*u+sin(thetas).*v;
m2 = -sin(thetas).*u+cos(thetas).*v;

f2 = figure;
f2.GraphicsSmoothing = 'on';
f2.Renderer = 'painters';
plot3(x,y,z,'linewidth',2,'color','black'); hold on;
quiv1 = quiver3(midpoint(1,:),midpoint(2,:),midpoint(3,:),...
    scale*m1(1,:),scale*m1(2,:),scale*m1(3,:),...
    'color','red','linewidth',1);
quiv2 = quiver3(midpoint(1,:),midpoint(2,:),midpoint(3,:),...
    scale*m2(1,:),scale*m2(2,:),scale*m2(3,:),...
    'color','blue','linewidth',1);
plot3(x,y,z,'linewidth',2,'color','black'); hold on;
hold off;
axis manual;
axis equal;
axis([-2 2 -2 2 0 4]);
set (gcf,'Position',[0,0,512,512]);

nSteps = 20000;
eps = 1e-6;
grad = zeros(1,n-1);
%thetas = mod(linspace(0,3*pi,n-1),2*pi);
beta = 1;
twistE = zeros(1,nSteps);
writerObj = VideoWriter('3c_material_frame.avi');
open(writerObj);
frame = getframe;
for iter=1:nSteps
    %%% YOUR CODE TO COMPUTE TWIST ENERGY AND GRADIENT %%%
    m = thetas(2:end)-thetas(1:end-1);
    m = mod(m,2*pi);
    m(m>pi) = m(m>pi)-pi;
    l = diff_len(1:end-1)+diff_len(2:end);
    twistE(iter) = 0.5*beta*sum(m.*m./l);
    quot = m./l;
    grad(2:end-1) = 2*beta*(quot(1:end-1)-quot(2:end));
    grad(end) = 2*beta*quot(end);
    grad(1) = -2*beta*quot(1);
    %%% END HOMEWORK PROBLEM %%%
    stepSize = 0.038;
    thetas = thetas-stepSize*grad;
    m1 = cos(thetas).*u+sin(thetas).*v;
    m2 = -sin(thetas).*u+cos(thetas).*v;
    quiv1.UData = scale*m1(1,:); quiv1.VData = scale*m1(2,:); quiv1.WData = scale*m1(3,:);
    quiv2.UData = scale*m2(1,:); quiv2.VData = scale*m2(2,:); quiv2.WData = scale*m2(3,:);
    %if iter == 1 || mod(iter,nSteps/4) == 0 
    %    saveas(gcf,['3c_grad_descent_iter-' int2str(iter) '.png']);
    %end
    drawnow;
    if mod(iter,10) == 0
        frame = getframe;
        writeVideo(writerObj,frame);
    end
end
close(writerObj);