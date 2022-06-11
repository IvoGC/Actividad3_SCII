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
Km = 7.53*10^-3;
Ts=1*10^-4; dtao=1*10^-5 ;KMAX=100000;T=Ts*KMAX;

%DEFINO MATRICES
%X=[ia ; tita ; w];
A=[-RA/LAA 0 -Km/LAA  ; 0 0 1 ; Ki_par/J 0 -Bm/J];
B=[1/LAA; 0; 0];
C=[0 1 0];
D=[0];

sys=ss(A,B,C,D);
sys_d=c2d(sys,Ts,'zoh');

Ad=sys_d.a;
Bd=sys_d.b;


Qd=diag([1/0.01 1 1/1000]) ;  Rd=1000;

[K,sol,e]=dlqr(Ad,Bd,Qd,Rd);
PLC=eig(Ad-Bd*K)

% %implementacion de funciones a usar
% tf=2; dt=1*10^-5; t=0:dt:(tf-dt); periodo=1;%[seg]
% torq=1.15*10^-3;
%
% Ref=pi/2*square(2*pi*t/periodo);%funcion de referencia que varia entre pi/2 y -pi/2
% TL=torq/2*square(2*pi*t/periodo)+torq/2;%Funcion torque que varia entre 0 y 1.15*10^-3
%X=[ia ; tita ; w];
TL=1.15*10^-3;
TL=0;

%Implementación del controlador en el modelo no lineal en tiempo continuo.
X=[0;1;0]; t=0;
ial=X(1); tital=X(2); wl=X(3);
ia=X(1); tita=X(2); w=X(3);h=Ts/20; u=[];i=1;
u_k(1)=0;

for ki=1:KMAX
    %ul=-K*X+ref*G;
    ul=-K*X;
    X=Ad*X+Bd*ul;
    ial(ki)=X(1);
    tital(ki)=X(2);
    wl(ki)=X(3);
    u_kl(ki)=ul;
end
tl=(0:KMAX-1)*Ts;
x=[0;1;0];
for ki=1:KMAX
    %X=[ia ; tita ; w];
    ud(ki)=-K*x;
    %ul=-K*X+ref*G; como seria la G en tiempo discreto
    for kii=1:Ts/h
        u(i)=ud(ki);
        ia_p=-RA/LAA*ia(i)-Km/LAA*w(i)+1/LAA*u(i);  %ia_p
        tita_p= w(i);                                  %tita_p
        w_p=Ki_par/J*ia(i)-Bm/J*w(i)-1/J*TL;        %W_pq
        
        ia(i+1)=ia(i)+h*ia_p;
        tita(i+1)=tita(i)+h*tita_p;
        w(i+1)=w(i)+h*w_p;

        i=i+1;
    end
    x=[ia(i) ; tita(i) ; w(i)];
end
u(i)=ud(ki);t=0:h:T;

figure(1);
subplot(2,2,1);plot(tl,ial,'b',t,ia,'r');grid on; title('Corriente');hold on;
subplot(2,2,2);plot(tl,tital,'b',t,tita,'r');grid on;title('Ángulo');hold on;
subplot(2,2,[3:4]);plot(tl,wl,'b',t,w,'r');grid on;title('Velocidad Agular');hold on;
figure(2);
plot(tl,u_kl,'b',t,u,'k');grid on;title('Acción de control');xlabel('Tiempo en Seg.');hold on;


