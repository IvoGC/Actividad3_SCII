clc; clear ;
%{
---------------------------------------------------------------------------
Implementar un sistema en variables de estado que controle el ángulo del motor, para
consignas de pi/2 y –pi/2 cambiando cada 500 milisegundos y que el TL de 1,5 10-3
aparece sólo para pi/2, para –pi/2 es nulo. Hallar el valor de integración
Euler adecuado.
---------------------------------------------------------------------------
%}
%implementacion de funciones a usar
% PARAMETROS DE SIMULACION
tm=1e-4; tf=2;
%deltat=1*10^-5 ;
KMAX=tf/tm;
%Tt=tm*KMAX;
dt=tm/10; 
t=0:dt:tf;
periodo=1;%[seg]
torq=1.15*10^-3;
Ref=pi/2*square(2*pi*t/periodo);%funcion de referencia que varia entre pi/2 y -pi/2
TL=torq/2*square(2*pi*t/periodo)+torq/2;%Funcion torque que varia entre 0 y 1.15*10^-3
Vh=tm/dt;u=[];i=1;  

% figure
% subplot(2,1,1);plot(t,Ref,'r');grid on; title('Referencia');hold on;
% subplot(2,1,2);plot(t,TL,'g');grid on;title('Torque TL');hold on;


%DEFINO PARAMETROS
LAA = 366*10^-6;
J = 5*10^-9;
RA = 55.6;
Bm = 0;
Ki_par = 7.49*10^-3;
Km = 7.53*10^-3;

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

Ad=sys_d.a; 
Bd=sys_d.b;
Cd=sys_d.c;

Ao=Ad'; 
%Bo=C(2,:)';
Bo=Cd';
Co=Bd';

Aa=[Ad zeros(3,1);-C(2,:)*Ad 1];
Ba=[Bd ; -C(2,:)*Bd];
Ca=[1 0];

%CONTROLABILIDAD
M=[Ba Aa*Ba Aa^2*Ba Aa^3*Ba Aa^4*Ba];
rango=rank(M)

Mo=[Bo Ao*Bo Ao^2*Bo Ao^3*Bo];
rango_o=rank(Mo)

%CONTROLADOR K e Ki
%para el calculo del mismo se utiliza el metodo LQR para lo cual definimos
%Q=[ia ; w; tita; Ki];
%Q=diag([1 1/1000 5000 0.2]); R=2000;
% Q=diag([1/0.01 1/1000 1 2]); R=200;
Q=diag([1 1/1000 1 10]); R=1e2;

%Q=diag([1000 1 900 20000]); R=20;
Ka=dlqr(Aa,Ba,Q,R);
K_i=-Ka(4); K=Ka(1:3);

%Qo=diag([1 1 10]); Ro=[100 , 100 ; 1 , 1000];
Qo=1e0*(diag([1 1 1])); Ro=1e3*([1 , 0 ; 0 , 1]);
Ko=dlqr(Ao,Bo,Qo,Ro);
%Ko=place(Ao,Bo,[0.999999 0.99998 0.99997]);


Referencia=pi/2;
%ITERACION 
 
Ve=0; Xf=[0 0 0]; u_k(1)=0; 
Ref1(1)=0;
X=zeros(3,round(tf/dt));    %X=[ia ; w ; tita ];
X(:,1)=[0;0;0];     %condiciones iniciales de X
err=0;
x_int=[0 0 0];
%Xhat=zeros(3,round(tf/dt)); %Xhat=[iahat ; wrhat ; titahat ]
Xhat=[0;0;0];  %condiciones iniciales de Xhat
Up(1)=0;
for ki=1:KMAX

    u_k(ki)=-K*Xhat+K_i*Ve;
    %u_k(ki)=-K*Xf'+K_i*Ve;
    Y=C*x_int';
    
    zona_m=0.5;
    if(u_k(ki)<=zona_m)
       if((u_k(ki)>=-zona_m))
        u_k(ki)=0;
       else
           u_k(ki)=u_k(ki)+zona_m;
       end
    else
        u_k(ki)=u_k(ki)-zona_m;
    end
   sat=15;
    if (u_k(ki)>=sat  )
       u_k(ki)=sat;
    end
    if (u_k(ki)<=-sat)
       u_k(ki)=-sat;
    end
    
    for kii=1:Vh
        %X_a=[X(1,i), X(2,i), X(3,i)] ; %X=[ia ; w ; tita ];
        X_a=X(:,i)';
        u(i)=u_k(ki);
        
        %sistema original
        Xp_1=-RA/LAA*X_a(1)-(Km/LAA)*X_a(2)+(1/LAA)*u(i);   %ia_p
        Xp_2=(Ki_par/J)*X_a(1)-(Bm/J)*X_a(2) -(1/J)*TL(i);    %w_p 
        Xp_3= X_a(2);                                       %tita_p
        Xp_a=[Xp_1 , Xp_2 , Xp_3];
        Xf=X_a+ dt*Xp_a;
        X(:,i+1)=Xf;     
        i=i+1;
    end
    %observador
    Yhat=C*Xhat;
    err=Y-Yhat;
    Xhat=u(i-1)*Bd+Ko'*err+Ad*Xhat;
    Ve=Ve+Ref(i)-C(2,:)*x_int';
    x_int=X(:,i)';
end
u(i)=u_k(ki); 


%  figure(2);hold on;
% subplot(2,2,1);plot(t,X(1,:),'g');grid on; title('Corriente ia');hold on;
% subplot(2,2,2);plot(t,X(3,:),'r');grid on;title('Ángulo tita');hold on;
% plot(t,Ref,'k');
% subplot(2,2,3);plot(t,X(2,:),'c');grid on;title('velocidad angular');hold on;
% subplot(2,2,4);plot(t,u);grid on;title('Acción de control');xlabel('Tiempo en Seg.');hold on;

 figure(2);hold on;
subplot(2,2,[1:2]);plot(t,X(3,:),'r');grid on;title('Ángulo tita');hold on;
plot(t,Ref,'k');
subplot(2,2,[3:4]);plot(t,u);grid on;title('Acción de control');xlabel('Tiempo en Seg.');hold on;