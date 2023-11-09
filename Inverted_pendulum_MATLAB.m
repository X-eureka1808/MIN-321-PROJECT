m = 1; M = 5; l = 2; g = 10; d = 1;dt=0.05; w0=0; theta0=0.3;x0=0;v0=0;u=10;
A=[0 1 0 0; 0 0 m*g/M 0; 0 0 0 1; 0 0 (m+M)*g/(M*l) 0];
eig(A)
B=[0; 1/M;0;1/(M*l)];
rank(ctrb(A,B));
eigs=[-0.1; -0.2; -2.3; -3.4];
k=place(A,B,eigs)
eig(A-B*k);

for t=0:dt:10
    clf;
    X=[x0-0;v0-0;theta0-0;w0-0];
    u=-k*X;%for controllability
    %u=0; (no control system)
    atheta=(u*cos(theta0)-m*l*sin(theta0)*cos(theta0)*w0*w0+(M+m)*g*sin(theta0))/(M*l+m*l*sin(theta0)*sin(theta0));
    w1=w0+atheta*dt;
    theta1=theta0+w1*dt;
    ax=(m*l*cos(theta0)*atheta-m*l*sin(theta0)*w0*w0)/(M+m);
    v1=v0+ax*dt;
    x1=x0+v1*dt;
    x=[x0-0.5 x0-0.5 x0+0.5 x0+0.5];
    y=[0 0.5 0.5 0];
    axis equal
    xm=[x0-l*sin(theta0)-0.5  x0-l*sin(theta0)-0.5  x0-l*sin(theta0)+0.5  x0-l*sin(theta0)+0.5];
    ym=[l*cos(theta0)-0.5  l*cos(theta0)+0.5  l*cos(theta0)+0.5  l*cos(theta0)-0.5];
    xc=[x0 -100; x0-l*sin(theta0) 200];
    yc=[0.25 -0.25; l*cos(theta0) -0.25];
    fill(x,y,'cyan');
    hold on;
    rectangle('Position', [x0-l*sin(theta0)-0.25, l*cos(theta0)-0.25, 0.5, 0.5], 'Curvature', [1, 1],'FaceColor','b');
    fill(xc,yc,'y');
    rectangle('Position', [x0-0.5, -0.25, 0.25, 0.25], 'Curvature', [1, 1],'FaceColor','b');
    rectangle('Position', [x0+0.25, -0.25, 0.25, 0.25], 'Curvature', [1, 1],'FaceColor','b');
    %rectangle('Position', [x0-l, -l, 2*l, 2*l], 'Curvature', [1, 1]);
    axis([-5 5 -5 5]);
    drawnow;
    pause(0.1);
    hold on;
    x0=x1;
    v0=v1;
    w0=w1;
    theta0=theta1;
end

% A=[0 1 0 0; 0 0 m*g/M 0; 0 0 0 1; 0 0 (m+M)*g/M*l 0];
% eig(A);
% B=[0; 1/M;0;1/M*l];
% rank(ctrb(A,B));
% eigs=[-1.1; -1.2; -1.3; -1.4];
% k=place(A,B,eigs);
% eig(A-B*k)
% %xdot=(A-Bk)x
% %k=place(A,B,eigs)
