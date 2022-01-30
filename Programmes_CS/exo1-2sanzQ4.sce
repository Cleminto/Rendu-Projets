clear all; 
clf();

N = 200; 
L = 100; 
h = L/N; 
x = 0:h:L-h;
c1 = 6/5;
p1 = 0.2;
c2 = 1/6;
p2 = 0.6;

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

dt = 0.5; 
T = 100;
t = 0;
k = 0;

u = (3*(c1*cosh(sqrt(c1)*(x - p1*L)/2)^(-2))+3*(c2*cosh(sqrt(c2)*(x - p2*L)/2)^(-2)))';

mat1 = (2*Id + dt*Axxx);
w = u;

while t <= T 
    
    for i = 1:200
        Tw = mat1\(2*Id*u - 0.5*dt*Ax*(w^2));
        w = Tw; // Pt fixe de picard
        
    end
    u = 2*w - u;
    
    
    
    drawlater
    clf()
    plot(x,u')
    a=get("current_axes")
    a.data_bounds=[0,-0.1;L,4];
    drawnow
    
    if k == 0 then
        xs2png(0,'exo2-sanzQ4_0.png');
    end
    if k == 50 then
        xs2png(0,'exo2-sanzQ4_1.png');
    end
    if k == 100 then
        xs2png(0,'exo2-sanzQ4_2.png');
    end
    if k == 150 then
        xs2png(0,'exo2-sanzQ4_3.png');
    end
    if k == 200 then
        xs2png(0,'exo2-sanzQ4_4.png');
    end
    t = t + dt;
    k = k+1; 
end
