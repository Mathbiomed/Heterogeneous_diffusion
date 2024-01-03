function [sol,mass,D,G,vector_field_x,vector_field_y,fick_x,fick_y, drift_x,drift_y]=diffusion2d_v6_high_hetero(p, D0, L, peri, Obs, n, alpha, max_it, ns)

sol=cell(1,max_it/ns); 
vector_field_x=cell(1,max_it/ns); 
vector_field_y=cell(1,max_it/ns); 
fick_x=cell(1,max_it/ns); 
fick_y=cell(1,max_it/ns); 

drift_x=cell(1,max_it/ns); 
drift_y=cell(1,max_it/ns); 


x=linspace(-1.1*L,1.1*L,n); y=x; [xx,yy]=meshgrid(y,x);
G=sqrt(xx.^2+yy.^2)<=L; G=G.*(1-Obs);
rr=sqrt(xx.^2+yy.^2)-L/3;
D=(1.0*(rr.^2)./((rr.^2)+L^2/9)+0.1).*(sqrt(xx.^2+yy.^2)>=(L/3)).*(sqrt(xx.^2+yy.^2)<=(L));

D=0.2/1.1^2*D.^2; % ABM to PDE D= stepsize^2;

p=G.*p; p=p/sum(sum(p));
sol{1}=p;
Dp=D.*p;

vector_field_y{1}=(G([2:end,end],:).*G).*(Dp([2:end,end],:)-Dp)+(G.*G([1,1:end-1],:)).*(Dp-Dp([1,1:end-1],:));
vector_field_x{1}=(G(:,[2:end,end]).*G).*(Dp(:,[2:end,end])-Dp)+(G.*G(:,[1,1:end-1])).*(Dp-Dp(:,[1,1:end-1]));
  % G ensures flux at boundary to be zero

temp_fick_y=(G([2:end,end],:).*G).*(p([2:end,end],:)-p)+(G.*G([1,1:end-1],:)).*(p-p([1,1:end-1],:));
fick_y{1}=D.*temp_fick_y;
temp_fick_x=(G(:,[2:end,end]).*G).*(p(:,[2:end,end])-p)+(G.*G(:,[1,1:end-1])).*(p-p(:,[1,1:end-1]));
fick_x{1}=D.*temp_fick_x;

temp_drift_y=(G([2:end,end],:).*G).*(D([2:end,end],:)-D)+(G.*G([1,1:end-1],:)).*(D-D([1,1:end-1],:));
drift_y{1}=p.*temp_drift_y;
temp_drift_x=(G(:,[2:end,end]).*G).*(D(:,[2:end,end])-D)+(G.*G(:,[1,1:end-1])).*(D-D(:,[1,1:end-1]));
drift_x{1}=p.*temp_drift_x;





mass=zeros(1,max_it+1);
mass(1)=sum(sum(p));
for it=1:max_it
    Dp=D.*p;
    p=p+alpha*((G([2:end,end],:).*G).*(Dp([2:end,end],:)-Dp) ...
        -(G.*G([1,1:end-1],: )).*(Dp-Dp([1,1:end-1],:)) ...
        +(G(:,[2:end,end]).*G).*(Dp(:,[2:end,end])-Dp) ...
        -(G.*G(:,[1,1:end-1])).*(Dp-Dp(:,[1,1:end-1])));
%     qy=0.5*(G([2:end,end],:)+G).*(Dp([2:end,end],:)-Dp)+0.5*(G+G([1,1:end-1],:)).*(Dp-Dp([1,1:end-1],:));
%     qx=0.5*(G(:,[2:end,end])+G).*(Dp(:,[2:end,end])-Dp)+0.5*(G+G(:,[1,1:end-1])).*(Dp-Dp(:,[1,1:end-1]));
    qy=(G([2:end,end],:).*G).*(Dp([2:end,end],:)-Dp)+(G.*G([1,1:end-1],:)).*(Dp-Dp([1,1:end-1],:));
    qx=(G(:,[2:end,end]).*G).*(Dp(:,[2:end,end])-Dp)+(G.*G(:,[1,1:end-1])).*(Dp-Dp(:,[1,1:end-1]));
    qy=qy/2; qx=qx/2;
    
    
    fy=(G([2:end,end],:).*G).*(p([2:end,end],:)-p)+(G.*G([1,1:end-1],:)).*(p-p([1,1:end-1],:));
    fx=(G(:,[2:end,end]).*G).*(p(:,[2:end,end])-p)+(G.*G(:,[1,1:end-1])).*(p-p(:,[1,1:end-1]));
    fy=D.*fy;
    fx=D.*fx;
    fy=fy/2; fx=fx/2;

    dy=(G([2:end,end],:).*G).*(D([2:end,end],:)-D)+(G.*G([1,1:end-1],:)).*(D-D([1,1:end-1],:));
    dx=(G(:,[2:end,end]).*G).*(D(:,[2:end,end])-D)+(G.*G(:,[1,1:end-1])).*(D-D(:,[1,1:end-1]));
    dy=p.*dy;
    dx=p.*dx;
    dy=dy/2; dx=dx/2;

    
    mass(it+1)=sum(sum(p));
    
    if mod(it,ns)==0
        sol{1+it/ns}=p;
        vector_field_x{1+it/ns}=qx;
        vector_field_y{1+it/ns}=qy;
        fick_x{1+it/ns}=fx;
        fick_y{1+it/ns}=fy;
        drift_x{1+it/ns}=dx;
        drift_y{1+it/ns}=dy;
    end
end
end