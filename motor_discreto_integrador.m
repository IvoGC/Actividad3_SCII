clc; clear ; close all
%{
---------------------------------------------------------------------------
Implementar un sistema en variables de estado que controle el ángulo del motor, para
consignas de pi/2 y –pi/2 cambiando cada 500 milisegundos y que el TL de 1,5 10-3
aparece sólo para pi/2, para –pi/2 es nulo. Hallar el valor de integración
Euler adecuado.
---------------------------------------------------------------------------
%}

%DEFINO PARAMETROS
LAA = 366*10^-6;
J = 5*10^-9;
RA = 55.6;
Bm = 0;
Ki_par = 7.49*10^-3;
Km = 7.53*10^-4;
% PARAMETROS DE SIMULACION
tm=1*10^-4; ts=1;
%deltat=1*10^-5 ;
dt=tm/20;
KMAX=ts/tm;
%Tt=tm*KMAX;

%DEFINO MATRICES
%X=[ia ; tita ; w];
A=[-RA/LAA 0 -Km/LAA  ; 0 0 1 ; Ki_par/J 0 -Bm/J];
B=[1/LAA; 0; 0];
C=[0 1 0];
D=[0];
I=eye(3);

sys=ss(A,B,C,D);
sys_d=c2d(sys,tm,'zoh');

Ad=sys_d.a;
Bd=sys_d.b;

Aa=[Ad zeros(3,1); -C*Ad 1];
Ba=[Bd ; -C*Bd];
Ca=[1 0];
%CONTROLABILIDAD
M=[Ba Aa*Ba Aa^2*Ba Aa^3*Ba Aa^4*Ba];
rango=rank(M)

%CONTROLADOR K e Ki
%para el calculo del mismo se utiliza el metodo LQR para lo cual definimos
%Q=[ia ; tita ; w; Ki];
Q=diag([1 1 1 1]); R=1;

Ka=dlqr(Aa,Ba,Q,R);
K_i=-Ka(4); K=Ka(1:3);

%ITERACION X=[ia ; tita ; w];
X=zeros(3,(tm/dt));
X(:,1)=[0,0,0]; ref=1;
ia=X(1); tita=X(2) ; wr=X(3); Vh=tm/dt;h=dt;u=[];i=1;u_k(1)=0;
Ve(1)=0;
for ki=1:KMAX
    Ve(ki+1)=Ve(ki)+ref-C*X(:,ki);
    u_k(ki)=-K*X(:,ki)+K_i*Ve(ki);
    Y=C*X(:,ki);
    for kii=1:Vh
        X_a=[X(1,i), X(2,i), X(3,i)] ; %X=[ia ; tita ; w];
        u(i)=u_k(ki);
        Xp_1=-RA/LAA*X_a(1)-Km/LAA*X_a(3)+1/LAA*u(i);   %ia_p
        Xp_2= X_a(3);                                   %tita_p
        Xp_3=Ki_par/J*X_a(1)-Bm/J*X_a(3)-1/J;           %W_p (TL(i))lo ignoramos momentaneamente
    
        Xp_a=[Xp_1 , Xp_2 , Xp_3];
    
        Xf=X_a+ dt*Xp_a;
        X(:,i+1)=Xf;
        i=i+1;
        
    end
    
    
end
u(i)=u_k(ki); t=0:h:KMAX*Vh*h;


figure(1);hold on;
subplot(2,2,1);plot(t,X(1,:),'r');grid on; title('Corriente ia');hold on;
subplot(2,2,2);plot(t,X(2,:),'g');grid on;title('Ángulo tita');hold on;
subplot(2,2,3);plot(t,X(3,:),'c');grid on;title('velocidad angular');hold on;
subplot(2,2,4);plot(t,u,'k');grid on;title('Acción de control');xlabel('Tiempo en Seg.');hold on;



