clc; clear ; close all
%{
---------------------------------------------------------------------------
 Calcular sistema controlador que haga evolucionar al péndulo en el equilibrio estable.
Objetivo de control: partiendo de una condición inicial nula en el desplazamiento 
y el ángulo en pi, hacer que el carro se desplace a 10 metros evitando las 
oscilaciones de la masa m, considerando que es una grúa. Una vez que delta=10
modificar a m a un valor 10 veces mayor y volver al origen evitando oscilaciones.
---------------------------------------------------------------------------
%}
%DEFINO PARAMETROS
m = 0.1; F = 0.1; l = 0.6; g = 9.8; M = 0.5; 

% PARAMETROS DE SIMULACION
tm=1e-3; tf=2;
%deltat=1*10^-5 ;
dt=1e-4;
KMAX=tf/dt;
%Tt=tm*KMAX;

%DEFINO MATRICES
%X=[delta delta_p phi phi_p]
%Sistema linealizado en el punto Xo
X_o=[0 0 0 0];%equilibio inestable phi=0
A = [0,1,0,0 ; 0,(-F/M),(-m*g/M),0 ; 0,0,0,1 ; 0,(-F/(l*M)),(((M+m)*-g)/(l*M)),0 ];
B = [0 ; (1/M) ; 0 ; (-1/(l*M))];
C = [1,0,0,0; 0,0,1,0];
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

Aa=[Ad zeros(4,1);-C(2,:)*Ad 1];
Ba=[Bd ; -C(2,:)*Bd];
Ca=[1 0];

%CONTROLABILIDAD
M_cont=[Ba Aa*Ba Aa^2*Ba Aa^3*Ba Aa^4*Ba Aa^5*Ba];
rango=rank(M_cont)

Mo=[Bo Ao*Bo Ao^2*Bo Ao^3*Bo Ao^4*Bo];
rango_o=rank(Mo)


%CONTROLADOR K e Ki
%para el calculo del mismo se utiliza el metodo LQR para lo cual definimos
%Q=[ia ; w; tita; Ki];
%Q=diag([1 1/1000 5000 0.2]); R=2000;
Q=diag([1000 1/100 1 1 0.2]); R=2000;
Ka=dlqr(Aa,Ba,Q,R);
K_i=-Ka(5); K=Ka(1:4);

%Qo=diag([1 1 10]); Ro=[100 , 100 ; 1 , 1000];
Qo=diag([1 1 10 1]); Ro=[100 , 1 ; 1 , 100];
Ko=dlqr(Ao,Bo,Qo,Ro);


%implementacion de funciones a usar
phi_i=pi;        %angulo    phi     inicial   ???
referencia1=10;  %posicion delta de referencia 1
referencia2=0;   %posicion delta de referencia 2
%generacion de funcion a usar
n=tf/dt;
Ref=zeros(1,round(n));  masa=zeros(1,round(n));

for j=1:1:n
    if (j < (n-1)/2 )
        Ref(j)=referencia1;
        masa(j)=m;
    else (j >= (n-1)/2 )
        Ref(j)=referencia2;
        masa(j)=m*10;
    end
end

Vh=tm/dt;u=[];i=1;  
t=0:dt:KMAX*(Vh)*dt;; 

% figure
% subplot(2,2,1);plot(Ref);title('REF');
% subplot(2,2,2);plot(masa);title('masa');

%ITERACION 
phi_i=pi;
d_i=0;
Ve(1)=0; Xf=[0 0 0 0]; u_k(1)=0; Xhatf=[0;0;0;0]; Xhat_a=[0;0;0;0];

X=zeros(4,round(tm/dt));   %X=[delta delta_p phi phi_p]
X(:,1)=[d_i,0,phi_i,0];     %condiciones iniciales de X

Xhat=zeros(4,round(tm/dt)); %X=[deltahat deltahat_p phihat phihat_p]
Xhat(:,1)=[d_i,0,phi_i,0];  %condiciones iniciales de Xhat

for ki=1:KMAX
    
    Ve(ki+1)=Ve(ki)+Ref(ki)-C(2,:)*Xf';
    u_k(ki)=-K*Xf'+K_i*Ve(ki);
    %u_k(ki)=-K*Xhatf+K_i*Ve(ki);
    Y=C*Xf';
    Yhat=C*Xhatf;
    err=Y-Yhat;
    
    for kii=1:Vh
        %X=[delta delta_p phi phi_p]        
        X_a=X(:,i)';
        u(i)=u_k(ki);
        Xp_1= X_a(2);                                   %delta_p
        %Xp_2=-F*X_a(2)/M-masa(i)*g*(X_a(3)-pi)/M+u(i)/M;  %delta_pp
        Xp_2=-F*X_a(2)/M-m*g*(X_a(3)-pi)/M+u(i)/M;  %delta_pp
        Xp_3= X_a(4);                                   %phi_p
        %Xp_4= -F*X_a(2)/(M*l)-g*(M+masa(i))*(X_a(3)-pi)/(M*l)-u(i)/(M*l);  %phi_pp
        Xp_4= -F*X_a(2)/(M*l)-g*(M+m)*(X_a(3)-pi)/(M*l)-u(i)/(M*l);  %phi_pp

        Xp_a=[Xp_1 , Xp_2 , Xp_3, Xp_4];
        
        Xf=X_a+ dt*Xp_a;
        X(:,i+1)=Xf;
        %Xhat_a=Xhat(:,i)';
        Xhat_a=[Xhat(1,i) ; Xhat(2,i) ; Xhat(3,i); Xhat(4,i)];
        Xhat_p=u(i)*B+Ko'*err+A*Xhat_a;
        Xhatf=Xhat_a + dt*Xhat_p;
        Xhat(:,i+1)=[Xhatf(1) ; Xhatf(2); Xhatf(3) ; Xhatf(4)];         
        i=i+1;
    end
    %Xhat(:,i-1)=Xhatf; 
end
u(i)=u_k(ki); 



%X=[delta delta_p phi phi_p]  
figure(1);hold on;
subplot(3,2,1);plot(t,X(1,:),'r');grid on; title('Distancia delta');hold on;
subplot(3,2,2);plot(t,X(3,:),'g');grid on;title('Delta punto');hold on;
subplot(3,2,3);plot(t,X(2,:),'c');grid on;title('Angulo Phi en el equilirbio estable');hold on;
subplot(3,2,4);plot(t,X(4,:),'r');grid on; title('Phi punto');hold on;
subplot(3,2,[5:6]);plot(t,u,'k');grid on;title('Acción de control');xlabel('Tiempo en Seg.');hold on;


