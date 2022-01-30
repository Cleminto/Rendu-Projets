clear all; 
clf();

N = 400; 
L = 50; 
h = L/N; 
x = 0:h:L-h;
c2 = 0.5;
d2 = 0.2*L;



// -- Derivee deuxieme -- // 

Axx = 2*eye(N,N) - diag(ones(N-1,1),1) - diag(ones(N-1,1),-1) 
Axx(1,N) = -1; 
Axx(N,1) = -1;

Axx = sparse((1/(h^2))*Axx);

// -- Derivee premiere -- // 

Ax = zeros(N,N) + diag(ones(N-1,1),1) - diag(ones(N-1,1),-1);
Ax(1,N) = -1;
Ax(N,1) = 1;
Ax = sparse((1/(2*h))*Ax);

// --  Schema en temps -- //  
Id = speye(N,N);

dt = 0.01; 
T = 100;
t = 0;
k = 0;
nper = 2500

inCos = sqrt((1-c2)/(2*c2))*(x-d2)
u0 = (-6*(c2 - 1)*(cosh(inCos)^(-2)) - 3*(1-c2))';
w = (-6*(c2 - 1)*(cosh(inCos)^(-2)) - 3*(1-c2))';



mat1 = (Id  + Axx);
mat2 = (Id + dt*Ax + Axx);


consMass = zeros(1,T-1);
consNorm = zeros(1,T-1);
while t <= T+1
    
    u1 = sum(u0);
    u1norm = sum(h*(u0^2)) + sum(((u0(2:N) - u0(1:N-1))^2)/h); 
    
    
    for i=1:100
        Tw = mat2\(mat1*u0 - 0.5*dt*Ax*(w^2))
        w = w + 0.2*(Tw - w);
    end
    u = w;
    
    u2 = sum(u);
    u2norm = sum(h*(u^2)) + sum(((u(2:N) - u(1:N-1))^2)/h);
    
//    if modulo(k,nper) == 0
//    drawlater
//    clf()
//    
//    
//  
//    drawnow
//    
//    end
    
//    if k == 0 then
//        xs2png(0,'exo3_0.png');
//    end
//    if k == 2500 then
//        xs2png(0,'exo3_1.png');
//    end
//    if k == 5000 then
//        xs2png(0,'exo3_2.png');
//    end
//    if k == 7500 then
//        xs2png(0,'exo3_3.png');
//    end
    if k == 10000 then
         plot(x,u')
        a=get("current_axes")
        a.data_bounds=[0,-1.6;L,3];
        xs2png(0,'exo3_4.png');
    end
    
    consMass(k+1) = u2 - u1;
    consNorm(k+1) = u2norm - u1norm;
    
    disp(k)
    t = t + dt; 
    k = k+1;
    
    u0 = u;
    
    
end


figure(0)
clf()
plot(0:max(size(consNorm))-1,consNorm)
a=get("current_axes")
a.data_bounds=[0,-1;T,1];
xs2png(0, 'exo3_Norm.png');


figure(0)
clf()
plot(0:max(size(consMass))-1,consMass)
a=get("current_axes")
a.data_bounds=[0,-1;T,1];
xs2png(0, 'exo3_Mass.png');
