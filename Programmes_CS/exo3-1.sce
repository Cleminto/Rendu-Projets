clear;

exec("CG.sce");

L=100;
N = 500;
h=L/(N-1);
x=0:h:L;







nu = 0.1; 

T=20;
dt=0.0001;
t=0;
k = 0;

// ----Matrice identité 

Id = speye(N,N);

//---Matrice de la dérivée première 
M1=zeros(N,N)-diag(ones(N-1,1),-1)+diag(ones(N-1,1),1);
M1(1,N)=-1;
M1(N,1)=1;
M1=sparse(M1/(2*h));


//----Matrice de la dérivée seconde 
M= 2*diag(ones(N,1))-diag(ones(N-1,1),1)-diag(ones(N-1,1),-1);

M(1,N)=-1;
M(N,1)=-1;
M2=sparse(M/h^2);











MAT = [Id+M2+dt*nu*M2 dt*M1; dt*M1 Id+dt*nu*M2];




eta0 = cos(8*%pi*x/L);
u0 = zeros(1,N);


X0 = [u0,eta0]';


X = X0;


f = ones(N,1)




while t<= T+dt
    
    F = [dt*f + (Id+M2)*X0(1:N) - (dt/2)*M1*(X0(1:N).^2);X(N+1:2*N)-dt*M1*(X0(N+1:2*N).*X0(1:N))]
    

    X = MAT\F;
    
    
    
    
    
    drawlater
    
    clf()
    subplot(2,1,1)
     plot(x,X(1:N)')
    a = gca();
    a.data_bounds=[0 -1.5;100 1.5]
    subplot(2,1,2)
    plot(x,X(N+1:2*N)')
    a = gca();
    a.data_bounds=[0 -2;100 2]

    
    
    
    drawnow
    
    
    
        
    X0(1:N) = X(1:N)
    X0(N+1:2*N) = X(N+1:2*N)
    
    
    
    t = t+dt;
    k = k+1;
    disp(t)
    
end
