clc; clear ; close all
%{
---------------------------------------------------------------------------
Implementar un sistema en variables de estado que controle el ángulo del motor, para
consignas de pi/2 y –pi/2 cambiando cada 500 milisegundos y que el TL de 1,5 10-3
aparece sólo para pi/2, para –pi/2 es nulo. Hallar el valor de integración
Euler adecuado.
---------------------------------------------------------------------------
%}

clc;clear all;
%DEFINO PARAMETROS
LAA = 366*10^-6;
J = 5*10^-9;
RA = 55.6;
Bm = 0;
Ki_par = 7.49*10^-3;
Km = 7.53*10^-3;

Ts=0.01;KMAX=5000;TEuler=1e-4;
TamanioFuente=20;
%Condiciones iniciales
alfa(1)=.1; color='.r';colorc='r';
% alfa(1)=.7; color='.b';colorc='b';
alfa(1)=.8; color='.b';colorc='k';
ref=10;
%DEFINO MATRICES
%X=[ia ; tita ; w];
Ac=[-RA/LAA 0 -Km/LAA  ; 0 0 1 ; Ki_par/J 0 -Bm/J];
Bc=[1/LAA; 0; 0];
Cc=[0 1 0];
Dc=[0];
I=eye(3);
sys_c=ss(Ac,Bc,Cc,Dc);
sys_d=c2d(sys_c,Ts,'zoh');
Ad=sys_d.a;
Bd=sys_d.b;
Aa=[Ad,zeros(4,1);-Mat_C*Ad, 1];
Ba=[Bd;-Mat_C*Bd];
Ma=[Ba Aa*Ba Aa^2*Ba Aa^3*Ba];%Matriz Controlabilidad
rango=rank(Ma);

%Cálculo del controlador por asignación de polos
auto_val=eig(Aa);
c_ai=poly(auto_val);
Wa=[c_ai(5) c_ai(4) c_ai(3) c_ai(2) 1;c_ai(4) c_ai(3) c_ai(2) 1 0;c_ai(3) c_ai(2) 1 0
0;c_ai(2) 1 0 0 0;1 0 0 0 0];
T=Ma*Wa;
A_controlable=inv(T)*Aa*T %Verificación de que T esté bien
%Ubicación de los polos de lazo cerrado en mui:
mui(1)=0.997; mui(2)=0.998; mui(3)= conj(mui(2)); mui(4)=0.99; mui(5)=0.99;
alfa_i=poly(mui);
Ka=fliplr(alfa_i(2:6)-c_ai(2:6))*inv(T);
K=Ka(1:4);KI=-Ka(5);
abs(eig(Aa-Ba*Ka))
%____________________________OBSERVADOR______________
Mat_M_Obs=[Mat_C;(Mat_C*Ad);(Mat_C*Ad^2);(Mat_C*Ad^3)];
rank(Mat_M_Obs)
%Como el Rango es 4, se genera el sistema Dual y se vuelve a repetir la
%secuencia de cálculo vista
Mat_C_O=Bd';
Mat_A_O=Ad';
Mat_B_O=Mat_C';
Mat_M_D=Mat_M_Obs';%[Mat_Ba Mat_Aa*Mat_Ba Mat_Aa^2*Mat_Ba Mat_Aa^3*Mat_Ba];%Matriz 
Controlabilidad
c_aid=poly(Ad);
Wa=[c_aid(4) c_aid(3) c_aid(2) 1;c_aid(3) c_aid(2) 1 0;c_aid(2) 1 0 0;1 0 0 0];
Mat_T_D=Mat_M_D*Wa;
alfa_i_O=poly([.0999 .095 .094 .093]); %Polos deseados del Observador muy bueno para el caso 
lineal
alfa_i_O=poly([.3 .5 .94 .3]); %Polos deseados del Observador
Ko=(fliplr(alfa_i_O(2:5)-c_aid(2:5))*inv(Mat_T_D))';
abs(eig(Ad-Ko*Mat_C))
t=0; x=[0;0;alfa(1);0];
p(1)=x(1); p_p(1)=x(2); alfa(1)=x(3); omega(1)=x(4);ve(1)=0;
u_k(1)=0;xang=[0 ;0;0;0];
for ki=2:KMAX
 t=[t ki*Ts];
 ve(ki)=ve(ki-1)+ref-Mat_C*x;
 % u=-K*x+KI*ve(ki);
 u=-K*xang+KI*ve(ki);
 ys=Mat_C*x; %Acá DEBE medirse y.
 x=Ad*x+Bd*u;
 p(ki)=x(1);
 p_p(ki)=x(2);
 alfa(ki)=x(3);
 omega(ki)=x(4);
 u_k(ki)=u;
 xang=Ad*xang+Bd*u+Ko*(ys-Mat_C*xang);%Acá se usa y.
end
u=u_k;figure(1);hold on;
subplot(3,2,1);plot(t,omega,color);grid on; title('Velocidad ángulo','FontSize',TamanioFuente);hold on;
subplot(3,2,2);plot(t,alfa,color);grid on;title('Ángulo','FontSize',TamanioFuente);hold on;
subplot(3,2,3); plot(t,p,color);grid on;title('Posición carro','FontSize',TamanioFuente);hold on;
subplot(3,2,4);plot(t,p_p,color);grid on;title('Velocidad carro','FontSize',TamanioFuente);hold on;
subplot(3,1,3);plot(t,u,color);grid on;title('Acción de control','FontSize',TamanioFuente);xlabel('Tiempo en Seg.','FontSize',TamanioFuente);hold on;
figure(2);hold on;
subplot(2,2,1);plot(alfa,omega,color);grid on;xlabel('Ángulo','FontSize',TamanioFuente);ylabel('Velocidad angular','FontSize',TamanioFuente);hold on;
subplot(2,2,2);plot(p,p_p,color);grid on;xlabel('Posición carro','FontSize',TamanioFuente);ylabel('Velocidad carro','FontSize',TamanioFuente);hold on;




%Verificación de la solución con el modelo no lineal en tiempo continuo.
T=t(ki);x=[0;0;alfa(1);0];
p=x(1); p_p=x(2); alfa=x(3); omega=x(4); tita_pp(1)=0;Vh=Ts/TEuler;h=TEuler; u=[];i=1;
u_k(1)=0;xang=[0;0;0;0];
for ki=1:KMAX
    ve(ki+1)=ve(ki)+ref-Mat_C*x;
    u1(ki)=-K*xang+KI*ve(ki); %Control con observación de estados
    % u1(ki)=-K*x+KI*ve(ki); %Sin Observador
    ys=Mat_C*x;%Acá se mide la salida.
    for kii=1:Vh
        u(i)=u1(ki);
        p_pp=(1/(M+m))*(u(i)-m*long*tita_pp*cos(alfa(i))+m*long*omega(i)^2*sin(alfa(i))-Fricc*p_p(i));
        tita_pp=(1/long)*(g*sin(alfa(i))-p_pp*cos(alfa(i)));
        p_p(i+1)=p_p(i)+h*p_pp;
        p(i+1)=p(i)+h*p_p(i);
        omega(i+1)=omega(i)+h*tita_pp;
        alfa(i+1)=alfa(i)+h*omega(i);
        i=i+1;
    end
    x=[p(i-1); p_p(i-1); alfa(i-1); omega(i-1)];
    xang=Ad*xang+Bd*u1(ki)+Ko*(ys-Mat_C*xang);%Acá se usa y.
end
u(i)=u1(ki);t=0:h:KMAX*(Vh)*h;
figure(1);hold on;
subplot(3,2,1);plot(t,omega,colorc);grid on; title('Velocidad ángulo','FontSize',TamanioFuente);hold on;
subplot(3,2,2);plot(t,alfa,colorc);grid on;title('Ángulo','FontSize',TamanioFuente);hold on;
subplot(3,2,3);plot(t,p,colorc);grid on;title('Posición carro','FontSize',TamanioFuente);hold on;
subplot(3,2,4);plot(t,p_p,colorc);grid on;title('Velocidad carro','FontSize',TamanioFuente);hold on;
subplot(3,1,3);plot(t,u,colorc);grid on;title('Acción de control','FontSize',TamanioFuente);xlabel('Tiempo en Seg.','FontSize',TamanioFuente);hold on;
figure(2);hold on;
subplot(2,2,1);plot(alfa,omega,colorc);grid on;xlabel('Ángulo','FontSize',TamanioFuente);ylabel('Velocidad angular','FontSize',TamanioFuente);hold on;
subplot(2,2,2);plot(p,p_p,colorc);grid on;xlabel('Posición carro','FontSize',TamanioFuente);ylabel('Velocidad carro','FontSize',TamanioFuente);hold on;




