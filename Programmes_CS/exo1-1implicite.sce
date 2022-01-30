clear all; 
clf();

N = 200; 
L = 50; 
h = L/N; 
x = (0:h:L-h);

// -- Derivee troisieme -- // 

Axxx = -2*diag(ones(N-1,1),1) + 2*diag(ones(N-1,1),-1) + diag(ones(N-2,1),2) - diag(ones(N-2,1),-2)
Axxx(1,N) = 2; Axxx(1,N-1) = -1; Axxx(2,N) = -1;
Axxx(N,1) = -2; Axxx(N,2) = 1; Axxx(N-1,1) = 1;

Axxx = sparse((1/(2*h^3))*Axxx);

// --  Schema en temps -- //  


Id = speye(N,N);
dt = 1; 
T = 100;
t = 0;
k = 0;
u = cos(8*%pi*x/L)';



M = (Id + dt*Axxx);
while t <= T 
    
    u = M\u;
    
    
    drawlater
    clf()
    plot(x,u')
    a=get("current_axes")
    a.data_bounds=[0,-1;L,1];
    drawnow
    if k == 0 then
        xs2png(0,'exo1-impl_0.png');
    end
    if k == 50 then
        xs2png(0,'exo1-impl_1.png');
    end
    if k == 100 then
        xs2png(0,'exo1-impl_2.png');
    end
    
    t = t + dt;
    k = k+1; 
end
