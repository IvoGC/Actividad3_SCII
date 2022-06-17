clc; clear ; close all
%{
---------------------------------------------------------------------------
Ítem [2] Para el caso del avión, emplear un tiempo de integración por Euler
adecuado y un tiempo de simulación de 70seg. Los parámetros son a=0.07; w=9;
b=5; c=150, hallar un controlador para que los polos de lazo cerrado se ubican 
en ui=-15?15j; -0.5?0.5j, para referencias de 100 y -100 metros en altura, 
ambas con alturas iniciales de -500 y 500.
---------------------------------------------------------------------------
%}
%implementacion de funciones a usar
% PARAMETROS DE SIMULACION
tm=1e-3; tf=70;
KMAX=tf/tm;
dt=tm/10; 
t=0:dt:tf;
Referencia=+100;
H_inicial=+500;
Ref=Referencia*ones(1,round(tf/dt)+1);
Vh=tm/dt;u=[];i=1;

%DEFINO PARAMETROS
w=9; a=0.07; b=5; c=150;

%DEFINO MATRICES
%X=[alfa phi phi_p h]
A = [ -a,a,0,0 ; 0,0,1,0 ; w^2,-w^2,0,0 ; c,0,0,0 ];
B = [ 0; 0 ; w^2*b ; 0 ];
C = [ 0,0,0,1 ; 0,1,0,0];
D = 0;

sys=ss(A,B,C,D);
sys_d=c2d(sys,tm,'zoh');

Ad=sys_d.a; 
Bd=sys_d.b;
Cd=sys_d.c;

Ao=Ad';
Bo=Cd';
Co=Bd';

Aa=[Ad zeros(4,1);-C(1,:)*Ad 1];
Ba=[Bd ; -C(1,:)*Bd];
Ca=[1 0];
%CONTROLABILIDAD
M=[Ba Aa*Ba Aa^2*Ba Aa^3*Ba Aa^4*Ba Aa^5*Ba];
rango=rank(M)

Mo=[Bo Ao*Bo Ao^2*Bo Ao^3*Bo Ao^4*Bo];
rango_o=rank(Mo)

%CONTROLADOR K e Ki
% en este caso vamos a utilizar LQR
%[alfa phi phi_p h K_i]
%Q=1e-3*(diag([1 10 1 1000 0.03])); R=2e-3;
Q=1e-15*diag([1 1 1/10 1/100000 500]);    R=1;

Ka=dlqr(Aa,Ba,Q,R);
K_i=-Ka(5); K=Ka(1:4);
eig(Aa-Ba*Ka)
%CONTORLADOR Ko

%Qo=1e0*(diag([1 1 1 1])); Ro=1e0*([1 , 0 ; 0 , 1]);
Qo=1*diag([1 1 1 1]);    Ro=[1 0; 0 1];
Ko=dlqr(Ao,Bo,Qo,Ro);

%ITERACION 
Ve=0; u_k(1)=0; 
X=zeros(4,round(tf/dt));    %X=[alfa phi phi_p h]
X(:,1)=[0;0;0;H_inicial];
err=0; Xhat=[0;0;0;0];
x_int=[0 0 0 0]; 

for ki=1:KMAX
    %u_k(ki)=-K*Xhat+K_i*Ve;
    u_k(ki)=-K*x_int'+K_i*Ve;
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
    for kii=1:Vh
       X_a=X(:,i)';
       u(i)=u_k(ki);
       %sistema original
       Xp_1=a*(X_a(2)-X_a(1));              %alfa_p
       Xp_2= X_a(3);                        %phi_p
       Xp_3=-w^2*(X_a(2)-X_a(1)-b*u(i));    %phi_pp
       Xp_4=c*X_a(1);                       %h_p
       Xp_a=[Xp_1 , Xp_2 , Xp_3 , Xp_4];
       X(:,i+1)=X_a+ dt*Xp_a;
       
       i=i+1;
    end
    %observador
    Yhat=C*Xhat;
    err=Y-Yhat;
    Xhat=u(i-1)*Bd+Ko'*err+Ad*Xhat;
    Ve=Ve+Ref(i)-C(1,:)*x_int';
    x_int=X(:,i)';
end
u(i)=u_k(ki); 



figure(2);hold on;
subplot(3,2,1);plot(t,X(1,:),'g');grid on; title('Ángulo alfa');hold on;
subplot(3,2,2);plot(t,X(2,:),'r');grid on;title('Ángulo phi');hold on;
subplot(3,2,3);plot(t,X(3,:),'c');grid on;title('Phi punto');hold on;
subplot(3,2,4);plot(t,X(4,:),'r');grid on;title('Altura H');hold on;
plot(t,Ref,'k');
subplot(3,2,[5:6]);plot(t,u);grid on;title('Acción de control');xlabel('Tiempo en Seg.');hold on;




