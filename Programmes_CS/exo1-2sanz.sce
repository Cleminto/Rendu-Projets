clear all; 
clf();

N = 1000; 
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

dt = 1; 
T = 100;
t = 0;
k = 0;
nper = 10;

u = cos(8*%pi*x/L)';

mat1 = (2*Id + dt*Axxx);
w = cos(8*%pi*x/L)';//zeros(N,1);

consMass = zeros(1,T-1);
consNorm = zeros(1,T-1);

while t <= T 
    u1 = sum(u);
    u1norm = sum(h*u^2)
    
    for i = 1:50
        Tw = mat1\(2*Id*u - 0.5*dt*Ax*(w^2));
        w = Tw; // Pt fixe de picard
        
    end
    u = 2*w - u;
    
    u2 = sum(u);
    u2norm = sum(h*u^2)
    
    consMass(k+1) = u2-u1;
    consNorm(k+1) = u2norm-u1norm;
    
    
    if modulo(k,nper) == 0
    drawlater
    clf()
    plot(x,u')
    a=get("current_axes")
    a.data_bounds=[0,-1;L,2.5];
    drawnow
    end
    
    if k == 0 then
        xs2png(0,'exo2-sanz_0.png');
    end
    if k == 25 then
        xs2png(0,'exo2-sanz_1.png');
    end
    if k == 50 then
        xs2png(0,'exo2-sanz_2.png');
    end
    if k == 75 then
        xs2png(0,'exo2-sanz_3.png');
    end
    if k == 100 then
        xs2png(0,'exo2-sanz_4.png');
    end
    t = t + dt; 
    k = k+1;
    disp(k)
end

figure(0)
clf()
plot(0:dt:T ,consMass);
a=get("current_axes")
a.data_bounds=[0,-1;T,1];
xs2png(0, 'exo2-sanz_mass.png');

figure(0)
clf()
plot(0:dt:T ,consNorm);
a=get("current_axes")
a.data_bounds=[0,-1;T,1];
xs2png(0, 'exo2-sanz_norm.png');

