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

dt = 0.01; 
T = 100;
t = 0;
k = 0;

u = cos(8*%pi*x/L)';

mat1 = (Id + 0.5*dt*Axxx);
mat2 = (Id - 0.5*dt*Axxx);

while t <= T+dt
    u = mat1\(mat2*u - 0.5*dt*Ax*(u^2));
    
    
    
    drawlater
    clf()
    plot(x,u')
    a=get("current_axes")
    a.data_bounds=[0,-1;L,2.5];
    drawnow
    if k == 0 then
        xs2png(0,'exo2-expl_0.png');
    end
    if k == 2500 then
        xs2png(0,'exo2-expl_1.png');
    end
    if k == 5000 then
        xs2png(0,'exo2-expl_2.png');
    end
    if k == 7500 then
        xs2png(0,'exo2-expl_3.png');
    end
    if k == 10000 then
        xs2png(0,'exo2-expl_4.png');
    end
    t = t + dt; 
    k = k+1;
    disp(k)
    disp(t)
    
end
