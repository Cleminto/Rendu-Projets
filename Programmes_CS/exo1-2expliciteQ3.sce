clear all; 
clf();

N = 200; 
L = 100; 
h = L/N; 
x = 0:h:L-h;

// -- Derivee troisieme -- // 

Axxx = -2*diag(ones(N-1,1),1) + 2*diag(ones(N-1,1),-1) + diag(ones(N-2,1),2) - diag(ones(N-2,1),-2)
Axxx(1,N) = 2; Axxx(1,N-1) = -1; Axxx(2,N) = -1;
Axxx(N,1) = -2; Axxx(N,2) = 1; Axxx(N-1,1) = 1;

Axxx = sparse((1/(2*h^3))*Axxx);

// -- Derivee premiere -- // 

Ax = diag(ones(N-1,1),1) - diag(ones(N-1,1),-1);
Ax(1,N) = -1;
Ax(N,1) = 1;

Ax = sparse((1/(2*h))*Ax);

// --  Schema en temps -- //  
Id = speye(N,N);

dt = 0.04; 
T = 100;
t = 0;
c1 = 1/3;
p1 = 0.4;
u = 3*(c1*cosh(sqrt(c1)*(x - p1*L)/2)^(-2))';

mat1 = (Id + 0.5*dt*Axxx);
mat2 = (Id - 0.5*dt*Axxx);

while t < T 
    u = mat1\(mat2*u - 0.5*dt*Ax*(u^2));
    
    
    
    drawlater
    clf()
    plot(x,u')
    a=get("current_axes")
    a.data_bounds=[0,-1;L,2.5];
    drawnow
    t = t + dt; 
end
