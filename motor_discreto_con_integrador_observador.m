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
% PARAMETROS DE SIMULACION
tm=1e-3; tf=2;
%deltat=1*10^-5 ;
dt=tm/10;
KMAX=tf/tm;
%Tt=tm*KMAX;

%DEFINO MATRICES
%X=[ia ; w ; tita ];
A=[-RA/LAA -Km/LAA  0  ; Ki_par/J -Bm/J 0; 0 1 0 ];
B=[1/LAA; 0; 0];
C=[0 1 0; 0 0 1];
%C=[ 0 0 1];
D=[0];
%Dejamos explicitas las matrices del sistema dual que es el "observado"
sys=ss(A,B,C,D);
sys_d=c2d(sys,tm,'zoh');

Ad=sys_d.a; Bd=sys_d.b;

Ao=Ad'; 
Bo=C(2,:)';  %Bo es la matriz de entrada del Dual, y tiene dos columnas.
Bo=C';  
Co=Bd';

Aa=[Ad zeros(3,1);-C(2,:)*Ad 1];
Ba=[Bd ; -C(2,:)*Bd]; %Acá es correcto porque controlas una sola salida
Ca=[1 0];

%CONTROLABILIDAD
M=[Ba Aa*Ba Aa^2*Ba Aa^3*Ba Aa^4*Ba];
rango=rank(M)

Mo=[Bo Ao*Bo Ao^2*Bo Ao^3*Bo];
rango_o=rank(Mo)

%CONTROLADOR K e Ki
%para el calculo del mismo se utiliza el metodo LQR para lo cual definimos
%Q=[ia ; w; tita; Ki];
Q=diag([1 1/1000 100 0.2]); R=2000;
Ka=dlqr(Aa,Ba,Q,R);
K_i=-Ka(4); K=Ka(1:3);

Qo=diag([1 1 1/1000]); Ro=10;
Ko=dlqr(Ao,Bo,Qo,Ro);


%implementacion de funciones a usar
Vh=tm/dt;u=[];i=1;  
dt=tm/100;; t=0:dt:KMAX*(Vh)*dt;; periodo=1;%[seg]
torq=1.15*10^-3;
Ref=pi/2*square(2*pi*t/periodo);%funcion de referencia que varia entre pi/2 y -pi/2
TL=torq/2*square(2*pi*t/periodo)+torq/2;%Funcion torque que varia entre 0 y 1.15*10^-3


%ITERACION 
 
Ve(1)=0; Xf=[0 0 0]; u_k(1)=0; Xhatf=[0;0 ;0]; Xhat_a=[0;0;0];

X=zeros(3,round(tm/dt));    %X=[ia ; w ; tita ];
X(:,1)=[0,0,0];     %condiciones iniciales de X
ia=X(1); tita=X(3) ; wr=X(2);
Xhat=zeros(3,round(tm/dt)); %Xhat=[iahat ; wrhat ; titahat ]
Xhat(:,1)=[0,0,0];  %condiciones iniciales de Xhat
iahat=X(1); titahat=X(3) ; wrhat=X(2);
for ki=1:KMAX
    %t=[t ki*tm];
    Ve(ki+1)=Ve(ki)+Ref(ki)-C(2,:)*Xf';
    %u_k(ki)=-K*Xf'+K_i*Ve(ki);
    u_k(ki)=-K*Xhatf+K_i*Ve(ki);
    Y=C*Xf';
    Yhat=C*Xhatf;
    err=Y-Yhat;
    
    for kii=1:Vh
        %X_a=[X(1,i), X(2,i), X(3,i)] ; %X=[ia ; w ; tita ];
        
        X_a=X(:,i)';
        u(i)=u_k(ki);
        Xp_1=-RA/LAA*X_a(1)-(Km/LAA)*X_a(2)+(1/LAA)*u(i);   %ia_p
        Xp_2=(Ki_par/J)*X_a(1)-(Bm/J)*X_a(2)-1/J*TL(i);               %w_p (-1/J*TL(i))lo ignoramos momentaneamente
        Xp_3= X_a(2);                                       %tita_p

        Xp_a=[Xp_1 , Xp_2 , Xp_3];
        
        Xf=X_a+ dt*Xp_a;
        X(:,i+1)=Xf;
        %Xhat_a=Xhat(:,i)';
        %Xhat_a=[Xhat(1,i) ; Xhat(2,i) ; Xhat(3,i)];
        Xhat_p=u(i)*B+Ko'*err(1)+Ko'*err(2)+A*Xhat(:,i);
        Xhatf=Xhat(:,i) + dt*Xhat_p
                  
        i=i+1;
    end
    Xhat(:,i-1)=Xhatf; 
end
u(i)=u_k(ki); 




figure(1);hold on;
subplot(2,2,1);plot(t,X(1,:),'r');grid on; title('Corriente ia');hold on;
subplot(2,2,2);plot(t,X(3,:),'g');grid on;title('Ángulo tita');hold on;
subplot(2,2,3);plot(t,X(2,:),'c');grid on;title('velocidad angular');hold on;
subplot(2,2,4);plot(t,u,'k');grid on;title('Acción de control');xlabel('Tiempo en Seg.');hold on;

