clear; 
clf();

N = 100; 
L = 50; 
h = L/N; 
x = 0:h:L-h;

// -- Derivee troisieme -- // 

Axxx = -2*diag(ones(N-1,1),1) + 2*diag(ones(N-1,1),-1) + diag(ones(N-2,1),2) - diag(ones(N-2,1),-2)
Axxx(1,N) = 2; Axxx(1,N-1) = -1; Axxx(2,N) = -1;
Axxx(N,1) = -2; Axxx(N,2) = 1; Axxx(N-1,1) = 1;

Axxx = sparse((1/(2*h^3))*Axxx);

// --  Schema en temps -- //  

dt = 0.01; 
T = 10;
t = 0;
k = 0
u = cos(8*%pi*x/L)';



while t < T 
    u = u - dt*Axxx*u;
   
    drawlater
    clf()
    plot(x,u')
    a=get("current_axes")
    a.data_bounds=[0,-1;L,1];
    drawnow
    if k == 0 then
        xs2png(0,'exo1-expl_0.png');
    end
    if k == 500 then
        xs2png(0,'exo1-expl_1.png');
    end
    if k == 1000 then
        xs2png(0,'exo1-expl_2.png');
    end
    t = t + dt; 
    k = k+1;
    disp(t)
end
