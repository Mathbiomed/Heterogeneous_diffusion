clear; clc; close all;
load('cmap')

K=101; % number of grids
L=10; % length of domain
% ip=rand(K); ip=ip/sum(sum(ip)); % initial condition (random)
ip=ones(K); ip=ip/sum(sum(ip));
% ip=zeros(K); ip(51,51)=1;

D0=5; % diffusion coefficient
peri=0.2; % perinucleus width
max_it=250000; % maximum iteration
ns=max_it/100; % iteration number of print
alp=0.01; % ratio of dt/dx^2

x=linspace(-L,L,K); y=linspace(-L,L,K);
[xx,yy]=meshgrid(y,x);

Obs=xx.^2+yy.^2<(L/3)^2; % set obstacles

tic
[p,mass,D,G,vec_x, vec_y,fick_x,fick_y,drift_x,drift_y]=diffusion2d_high_hetero(ip, D0, L, peri, Obs, K, alp, max_it, ns);
toc
%% Heatmap for protein distribution over cytoplasm

fig1=figure;
set(gcf,'color','w');
imagesc(p{101})

colormap(cmap)


xticks([])
yticks([])
drawnow;   
caxis([0e-4,4e-4])


fig2=figure;
set(gcf,'color','w');
imagesc(p{1})
colormap(cmap)
%    shading interp; 
% axis([-L L -L L]);
xticks([])
yticks([])
drawnow;   
caxis([1.2e-4,2.5e-4])
     caxis([0e-4,4e-4])
  
%
fig3=figure;
set(gcf,'color','w');
imagesc(p{51})
colormap(cmap)
%    shading interp; 
% axis([-L L -L L]);
xticks([])
yticks([])
caxis([1.2e-4,2.5e-4])
caxis([0e-4,4e-4])

drawnow;   
       

% making short vectors (just for visualization)
% for i=1:101
%     vec_x= 0.25.* vec_x{i};
%     vec_y= 0.25.* vec_y{i};
% end

%% Drawing for flux for each time point

x_coarse=linspace(-L,L,11);
y_coarse=linspace(-L,L,11);
[xx_plot,yy_plot]=meshgrid(x_coarse, y_coarse);

scaling = 1*10^6; % Length scale 

figure
for i=1:5
    subplot(3,5,i)
    quiver(xx_plot, yy_plot, -vec_x{i*20-19}(1:10:101,1:10:101)*scaling,-vec_y{i*20-19}(1:10:101,1:10:101)*scaling,'AutoScale','off','LineWidth',1.5,'Color',[242,150,2]/255);
    xlim([-L,L])
    ylim([-L,L])
    xticks([])
    yticks([])
    
    
    subplot(3,5,i+5)
    quiver(xx_plot, yy_plot, -fick_x{i*20-19}(1:10:101,1:10:101)*scaling,-fick_y{i*20-19}(1:10:101,1:10:101)*scaling,'AutoScale','off','LineWidth',1.5,'Color',[141,192,30]/255);
    xlim([-L,L])
    ylim([-L,L])
    xticks([])
    yticks([])
    
    subplot(3,5,i+10)
    quiver(xx_plot, yy_plot, -drift_x{i*20-19}(1:10:101,1:10:101)*scaling,-drift_y{i*20-19}(1:10:101,1:10:101)*scaling,'AutoScale','off','LineWidth',1.5,'Color',[0,159,232]/255);
    xlim([-L,L])
    ylim([-L,L])
    xticks([])
    yticks([])
end