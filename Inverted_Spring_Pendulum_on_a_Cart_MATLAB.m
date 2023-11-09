m = 1;
M = 5; 
l = 3;
g = 10;
dt=0.01;
k=10;
w0=0; theta0=0.2; x0=0; v0=0; rdot0=0; r0=l;


A=[0 1 0 0; 0 0 m*g/M 0; 0 0 0 1; 0 0 ((M+m)*g*k)/(M*l*k-M*m*g) 0];
eig(A);
B=[0; 1/M; 0; 1/(M*(l-m*g/k))];
rank(ctrb(A,B));
eigs=[-1.1; -1.2; -1.3; -1.4];
d=place(A,B,eigs)
eig(A-B*d);

for t=0:dt:20
    clf;
    X=[x0+10;v0-0;theta0-0;w0-0];
    u=-d*X;
    ax=(u+k*(l-r0)*sin(theta0))/M;
    atheta=(u*cos(theta0)+k*(l-r0)*sin(theta0)*cos(theta0)+M*g*sin(theta0)-2*M*rdot0*w0)/(M*r0);
    ar=u*sin(theta0)/M+k*(l-r0)*sin(theta0)*sin(theta0)/M+w0*w0*r0-g*cos(theta0)+k*(l-r0)/m;
    w1=w0+atheta*dt;
    theta1=theta0+w1*dt;
    v1=v0+ax*dt;
    x1=x0+v1*dt;
    rdot1= rdot0+ar*dt;
    r1=r0+rdot1*dt;
    x2=[x0-0.5 x0-0.5 x0+0.5 x0+0.5];
    y2=[0 0.5 0.5 0];
    xm=[x0-r0*sin(theta0)-0.25 x0-r0*sin(theta0)-0.25 x0-r0*sin(theta0)+0.25 x0-r0*sin(theta0)+0.25];
    ym=[r0*cos(theta0)-0.25 r0*cos(theta0)+0.25 r0*cos(theta0)+0.25 r0*cos(theta0)-0.25];
    fill(x2,y2,[0.3010 0.7450 0.9330]);
    hold on;
    axis equal;
    rectangle('Position', [x0-r0*sin(theta0)-0.25, r0*cos(theta0)-0.25, 0.5, 0.5], 'Curvature', [1, 1],'FaceColor',"#A2142F");
    xc=[x0 -500; x0-r0*sin(theta0) 500];
    yc=[0.25 -0.25; r0*cos(theta0) -0.25];
    rectangle('Position', [x0-0.5, -0.25, 0.25, 0.25], 'Curvature', [1, 1],'FaceColor',"#0072BD");
    rectangle('Position', [x0+0.25, -0.25, 0.25, 0.25], 'Curvature', [1, 1],'FaceColor',"#0072BD");
    fill(xc,yc,'b');
    axis([-20 20 -10 10])

    drawnow;
    pause(0.1);
    x0=x1;
    v0=v1;
    r0=r1;
    theta0=theta1;
    w0=w1;
    rdot0=rdot1;
end

