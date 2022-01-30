clear;

exec("CG.sce");

L=100;
N = 120;
h=L/(N-1);
x=(-L/2):h:(L/2)+0.01;
y= x;

[X,Y] = meshgrid(x,y);




nu = 1; 

T=20;
dt=0.01;
t=0;
k = 0;

// ----Matrice identité 

Id = speye(N,N);

//---Matrice de la dérivée première 
M1=zeros(N,N)-diag(ones(N-1,1),-1)+diag(ones(N-1,1),1);
M1(1,N)=-1;
M1(N,1)=1;
M1=sparse(M1/(2*h));

Dx = Id.*.M1;
Dy = M1.*.Id;

DIV = Dx + Dy;


//----Matrice de la dérivée seconde 
M= 2*diag(ones(N,1))-diag(ones(N-1,1),1)-diag(ones(N-1,1),-1);

M(1,N)=-1;
M(N,1)=-1;
M = M/h^2;
M2=sparse(M.*.Id + Id.*.M);

nt = max(size(M2));

I = speye(nt,nt);

x2D = matrix(X,nt)
y2D = matrix(Y,nt)


b = 1;
d = 1;

MAT1 = (I + b*M2);
MAT2 = (I + d*M2);

A = 0.1;
Cs = 1+(A/2);

MAT = [MAT1 dt*Dx dt*Dy; dt*Dx MAT2 zeros(nt,nt); dt*Dy zeros(nt,nt) MAT2];


//eta0 = (A*(0.5*sqrt(3*A/Cs)*(y2D +10)));
//u0 = zeros(nt,1);
//v0 = eta0 - (eta0^2)/4;

eta0 = 0.5*exp(-(x2D.^2 + y2D.^2)/2);
u0 = zeros(nt,1);
v0 = zeros(nt,1);







eta = []
u = []
v = []




while t<= T+dt
    
    [eta,R] = CG(MAT1,MAT1*eta0 - dt*(Dx*u + Dy*v + Dx*(eta0.*u) + Dy*(eta0.*v)),eta0,1e-12,100); 
    [u,R] = CG(MAT2,MAT2*u0 - dt*(Dx*eta + 0.5*Dx*(u0.^2 + v0.^2)),u0,1e-12,100); 
    [v,R] = CG(MAT2,MAT2*v0 - dt*(Dy*eta + 0.5*Dy*(u0.^2 + v0.^2)),v0,1e-12,100); 
    
    //eta = MAT1\(MAT1*eta - dt*(Dx*u + Dy*v + Dx*(eta.*u) + Dy*(eta.*v)));
    //u = MAT2\(MAT2*u0 - dt*(Dx*eta + 0.5*Dx*(u0.^2 + v0.^2)));
    //v = MAT2\(MAT2*v0 - dt*(Dy*eta + 0.5*Dy*(u0.^2 + v0.^2)));
    
    
    
    eta2 = matrix(eta,N,N);
    drawlater
    clf()
    scf(0).color_map = jetcolormap(60);
    if max(eta2) > abs(min(eta2)) then
        Sgrayplot(x,y,eta2, rect = [-(L/2),-(L/2),L/2,L/2],zminmax = [-max(eta2),max(eta2)]);
        colorbar(-max(eta2),max(eta2))
    else 
        Sgrayplot(x,y,eta2, rect = [-(L/2),-(L/2),L/2,L/2],zminmax = [min(eta2),abs(min(eta2))]);
        colorbar(min(eta2),abs(min(eta2)))
    end
    
  
//    scf(0).color_map = jetcolormap(30)
//    surf(X,Y,eta2);
//    a = gca();
//    a.data_bounds=[-(L/2) -(L/2) -0.5;(L/2) (L/2) 0.5]
    
    
    
    drawnow
    
      if k == 0 then
       xs2png(gcf(),'exo4_2-0L100.png')
    elseif k == 500 then
       xs2png(gcf(),'exo4_2-5L100.png')
    elseif k == 1000 then
       xs2png(gcf(),'exo4_2-10L100.png')
    elseif k == 1500 then
       xs2png(gcf(),'exo4_2-15L100.png')
    elseif k == 2000 then
       xs2png(gcf(),'exo4_2-20L100.png')
    end
    
        
    eta0 = eta;
    u0 = u;
    v0 = v;
    
    
    t = t+dt;
    k = k+1;
    disp(t)
    
end
