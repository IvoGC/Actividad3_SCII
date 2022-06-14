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
%DEFINO PARAMETROS
w=3; a=0.05; b=5; c=100;

% PARAMETROS DE SIMULACION
tm=1e-2; tf=3.5;
%deltat=1*10^-5 ;
dt=tm/20;
KMAX=tf/dt;
%Tt=tm*KMAX;

%DEFINO MATRICES
%X=[alfa phi phi_p h]
A = [ -a,a,0,0 ; 0,0,1,0 ; w^2,-w^2,0,0 ; c,0,0,0 ];
B = [ 0; 0 ; w^2*b ; 0 ];
C = [ 0,0,0,1 ; 0,1,0,0];
D = 0;
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

Aa=[Ad zeros(4,1);-C(1,:)*Ad 1];
Ba=[Bd ; -C(1,:)*Bd];
Ca=[1 0];

%CONTROLABILIDAD
M=[Ba Aa*Ba Aa^2*Ba Aa^3*Ba Aa^4*Ba Aa^5*Ba];
rango=rank(M)

Mo=[Bo Ao*Bo Ao^2*Bo Ao^3*Bo Ao^4*Bo];
rango_o=rank(Mo)

%CONTROLADOR K, Ko e Ki 
%X=[alfa phi phi_p h]
%en este caso vamos a utilizar los plos brindados por la consigna
p1=-15+15i; p2=-15-15i; p3=-0.5+0.5i; p4=-0.5-0.5i;
%estos polos son del sistema original por lo que vamos a obtener los polos
%del sistema discreto referenciados a la circunferencia unitaria
pd1=exp(p1*tm); pd2=exp(p2*tm); pd3=exp(p3*tm); pd4=exp(p4*tm); 

%syms pd1 pd2 pd3 pd4 s
%expand((s-pd1)*(s-pd2)*(s-pd3)*(s-pd4))
a0=1;
a1=-p1-p2-p3-p4;
a2=p1*p2+p1*p3+p1*p4+p2*p3+p2*p4+p3*p4;
a3=-p1*p2*p3-p1*p2*p4-p1*p3*p4-p2*p3*p4;
a4=p1*p2*p3*p4;
phi_A=a0*A^4+a1*A^3+a2*A^2+a3*A^1+a4*A^0;

AUX=[B A*B A^2*B A^3*B];
K=[0 0 0 1]*inv(AUX)*phi_A;
K_i=-1/10000;
%preguntar sobre KI que pasa? donde va el polo

%ahora para el observador
po1=-50; po2=-30; po3=-1+i; po4=-1-i;
pod1=exp(po1*tm); pod2=exp(po2*tm); pod3=exp(po3*tm); pod4=exp(po4*tm); 

polos_o=[pod1 pod2 pod3 pod4];
Ko=place(Ao, Bo, polos_o);

%configuracion de parametros de simulacion
Vh=tm/dt;u=[];i=1;  
t=0:dt:KMAX*(Vh)*dt;; periodo=1;%[seg]
%CASO ALTURA INICIAL DE -100 Y ALTURA FINAL DE 500
H_i=100; H_f=-500; 
Ref=H_f;


%ITERACION 
 
Ve(1)=0; Xf=[0,0,0,0]; u_k(1)=0; Xhatf=[0;0;0;0]; Xhat_a=[0;0;0;0];

X=zeros(4,round(tm/dt));    %X=[ia ; w ; tita ];
X(:,1)=[0,0,0,H_i];     %condiciones iniciales de X
%OBSERVADOR
Xhat=zeros(4,round(tm/dt)); %Xhat=[iahat ; wrhat ; titahat ]
Xhat(:,1)=[0,0,0,H_i];  %condiciones iniciales de Xhat

for ki=1:KMAX
    Ve(ki+1)=Ve(ki)+Ref-C(1,:)*Xf';
    u_k(ki)=-K*Xf'+K_i*Ve(ki);
    %u_k(ki)=-K*Xhatf+K_i*Ve(ki);
    Y=C*Xf';
    Yhat=C*Xhatf;
    err=Y-Yhat;
    
    for kii=1:Vh
        %X=[alfa phi phi_p h]
        X_a=X(:,i)';
        u(i)=u_k(ki);
        Xp_1=a*(X_a(2)-X_a(1));             %alfa_p
        Xp_2=X_a(3);                         %phi_p
        Xp_3= -w^2*(X_a(2)-X_a(1)-b*u(i));  %phi_pp
        Xp_4= c*X_a(1);                     %h_p

        Xp_a=[Xp_1 , Xp_2 , Xp_3 ,Xp_4];
        
        Xf=X_a+ dt*Xp_a;
        X(:,i+1)=Xf;
        %Xhat_a=Xhat(:,i)';
        Xhat_a=[Xhat(1,i) ; Xhat(2,i) ; Xhat(3,i) ;Xhat(4,i)];
        Xhat_p=u(i)*B+Ko'*err+A*Xhat_a;
        Xhatf=Xhat_a + dt*Xhat_p;
        Xhat(:,i+1)=[Xhatf(1) ; Xhatf(2); Xhatf(3); Xhatf(4)];         
        i=i+1;
    end 
end
u(i)=u_k(ki); 




figure(1);hold on;
subplot(3,2,1);plot(t,X(1,:),'r');grid on; title('Angulo alfa');hold on;
subplot(3,2,2);plot(t,X(3,:),'g');grid on;title('Ángulo phi');hold on;
subplot(3,2,3);plot(t,X(2,:),'c');grid on;title('phi_punto');hold on;
subplot(3,2,4);plot(t,X(4,:),'r');grid on; title('Altura H');hold on;
subplot(3,2,[5:6]);plot(t,u,'k');grid on;title('Acción de control');xlabel('Tiempo en Seg.');hold on;

