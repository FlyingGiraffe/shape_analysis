clear; clc; close all;
%% Problem 2.c
a = 4; b = 2;
delta = pi/3;
n = 300;

t = linspace(0,2*pi,n);
x = sin(a*t*delta);
y = sin(b*t);

% u = zeros(1,n-2);
% v = zeros(1,n-2);
xy = [x; y];
diff = xy(:,2:end)-xy(:,1:end-1);

%%% YOUR CODE HERE TO COMPUTE GRADIENT HERE %%%
diff_len = sqrt(sum(diff.*diff));
diff = diff./diff_len; % normalize
uv = -diff(:,2:end)+diff(:,1:end-1);
u = uv(1,:);
v = uv(2,:);
%%% END HOMEWORK PROBLEM %%%

f1 = figure;
f1.GraphicsSmoothing = 'on';
f1.Renderer = 'painters';
plot(x,y,'linewidth',2,'color','black'); hold on;
quiver(x(2:end-1),y(2:end-1),u,v,10,'linewidth',1,'color','red');
axis equal;
saveas(gcf,'2c_derivative.png');

%% Problem 2.d
%%% YOUR CODE TO COMPUTE KAPPA HERE %%%
kappa = 4*sqrt(sum(uv.*uv))./(abs(diff_len(2:end))+abs(diff_len(1:end-1)));
%%% END HOMEWORK PROBLEM %%%

f2 = figure;
f2.GraphicsSmoothing = 'on';
f2.Renderer = 'painters';
X = x(2:end-1);
Y = y(2:end-1);
Z = zeros(size(kappa));
surface([X;X],[Y;Y],[Z;Z],[kappa;kappa],...
    'edgecolor', 'interp', 'linewidth',2);
axis equal;
colorbar;
saveas(gcf,'2d_curvature.png');

%% Problem 2.e
t0 = 0;
t1 = pi*1.25;
nSamples = 100;
nSteps = 5000;

% We provide a few examples of curves to try
%curveFunction = @(t) [(cos(t)-cos(3*t).^3); (sin(t)-sin(3*t).^3)]';
curveFunction = @(t) [cos(t);sin(t)]';
%curveFunction = @(t) [t;(t-t0).*(t1-t)]';
curve = curveFunction(linspace(t0,t1,nSamples));

colorFunction = @(i) [0.5-i/(nSteps*2);0.5-i/(nSteps*2);0.5-i/(nSteps*2)];

% Time step
%writerObj = VideoWriter('2d_grad_descent_1_consth.avi');
%open(writerObj);
f = figure;
plt = plot(curve(:,1),curve(:,2),'k','linewidth',2,'color',colorFunction(0)); hold on;
%plt = plot(curve(:,1),curve(:,2),'k','linewidth',2); hold on;
axis equal;
axis manual;
for i=1:nSteps
    %%% YOUR CODE HERE TO PERFORM GRADIENT DESCENT %%%
    d = curve(2:end,:)-curve(1:end-1,:);
    d = d./sqrt(sum(d.*d,2)); % normalize
    uv = zeros(nSamples,2);
    uv(2:end-1,:) = -d(2:end,:)+d(1:end-1,:);
    pow = 1;
    curve = curve-0.015*uv;
    %%% END HOMEWORK PROBLEM %%%
    if mod(i,nSteps/5) == 0
        plot(curve(:,1),curve(:,2),'k','linewidth',2,'color',colorFunction(i)); hold on;
    end
    %plt.XData = curve(:,1);
    %plt.YData = curve(:,2);
    %drawnow;
    %frame = getframe;
    %writeVideo(writerObj,frame);
end
saveas(gcf,['2e_grad_descent_1_nStep-' int2str(nSteps) '_consth.png']);
%saveas(gcf,['2e_grad_descent_3_nStep-' int2str(nSteps) '_pow-' num2str(pow,3) '.png']);
%close(writerObj);
