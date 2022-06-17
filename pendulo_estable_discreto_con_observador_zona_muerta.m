clc; clear ;
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
m = 0.1; F = 0.1; l = 2.6; g = 9.8; M = 0.5; 

% PARAMETROS DE SIMULACION
tm=1e-3; tf=15;
KMAX=tf/tm;
dt=tm/10; 
t=0:dt:tf;
Vh=tm/dt;u=[];i=1;
%FUNCIONES A USAR
referencia1=10;  %posicion delta de referencia 1
referencia2=0;   %posicion delta de referencia 2
%generacion de funcion a usar
n=round(tf/dt);Ref=zeros(1,n);  masa=zeros(1,n); 
for j=1:1:n+1
    if (j < (n)/2 )
        Ref(j)=referencia1;
        masa(j)=m;
    else if(j >= (n)/2 )
        Ref(j)=referencia2;
        masa(j)=m*10;
        end
    end
end

% figure
% subplot(2,1,1);plot(t,Ref,'g');title('Referencia');
% subplot(2,1,2);plot(t,masa,'r');title('Masa');

%DEFINO MATRICES
%X=[delta delta_p phi phi_p]
%Sistema linealizado en el punto Xo
X_o=[0 0 pi 0];%equilibio estable phi=pi
A = [0,1,0,0 ; 0,(-F/M),(-m*g/M),0 ; 0,0,0,1 ; 0,(-F/(l*M)),(((M+m)*-g)/(l*M)),0 ];
B = [0 ; (1/M) ; 0 ; (-1/(l*M))];
C = [1,0,0,0; 0,0,1,0];
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
M_cont=[Ba Aa*Ba Aa^2*Ba Aa^3*Ba Aa^4*Ba Aa^5*Ba];
rango=rank(M_cont)

Mo=[Bo Ao*Bo Ao^2*Bo Ao^3*Bo Ao^4*Bo];
rango_o=rank(Mo)

%CONTROLADOR K e Ki
% en este caso vamos a utilizar LQR
%[delta delta_p phi phi_p]
Q=1e0*(diag([1 10 150 1000 30])); R=2e3;
Ka=dlqr(Aa,Ba,Q,R);
K_i=-Ka(5); K=Ka(1:4);
%eig(Aa-Ba*Ka)

%CONTORLADOR Ko
Qo=1e0*(diag([1 1 1 1])); Ro=1e0*([1 , 0 ; 0 , 1]);
Ko=dlqr(Ao,Bo,Qo,Ro);
%eig(Ao-Bo*Ko)

%ITERACION 
Ve=0; u_k(1)=0; 
X=zeros(4,round(tf/dt));    %X=[delta delta_p phi phi_p]
X(:,1)=[0;0;0;0];           %condiciones iniciales
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
        X_a=X(:,i)'; %X=[delta delta_p phi phi_p]
        u(i)=u_k(ki);
        %sistema original
        Xp_1=X_a(2);
        Xp_2=-F*X_a(2)/M-masa(i)*g*(X_a(3)-pi)/M+u(i)/M;
        Xp_3=X_a(4);
        Xp_4=-F*X_a(2)/(M*l)-g*(M+masa(i))*(X_a(3)-pi)/(M*l)-u(i)/(M*l);
        Xp_a=[Xp_1 , Xp_2 , Xp_3, Xp_4];
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

% %X=[delta delta_p phi phi_p]
% figure(2);hold on;
% subplot(3,2,1);plot(t,X(1,:),'g');grid on; title('Distancia Delta');hold on;
% plot(t,Ref,'k');
% subplot(3,2,2);plot(t,X(2,:),'r');grid on;title('Delta Punto');hold on;
% subplot(3,2,3);plot(t,X(3,:),'c');grid on;title('Angulo phi');hold on;
% subplot(3,2,4);plot(t,X(4,:),'r');grid on;title('phi punto');hold on;
% subplot(3,2,[5:6]);plot(t,u);grid on;title('Acción de control');xlabel('Tiempo en Seg.');hold on;


 figure(2);hold on;
subplot(3,2,[1:2]);plot(t,X(1,:),'r');grid on; title('Distancia Delta');hold on;
plot(t,Ref,'k');
subplot(3,2,[3:4]);plot(t,X(3,:),'c');grid on;title('Angulo phi');hold on;
subplot(3,2,[5:6]);plot(t,u);grid on;title('Acción de control');xlabel('Tiempo en Seg.');hold on;


