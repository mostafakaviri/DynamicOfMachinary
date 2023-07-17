clc
close all
clear all
global l1 l2 l3 l4 l5 lBE h X m2 m3 m4 m5 m6  I2 I3 I4 I5 g
l1=21e-2; l2=15e-2; l3=18e-2; l4=20e-2; l5=25e-2; lBE=7.5e-2; h=24e-2; X=(25e-2)/3;
m2=5; m3=6; m4=7; m5=5; m6=2;
I2=1/12*m2*l2^2; I3=1/12*m3*l3^2; I4=1/12*m4*l4^2; I5=1/12*m5*l5^2;
g=9.81;

t_2=pi/4
dt_2=-8
ddt2=-5
i=0
time=.6
dt =.01
for t=0:dt:time
    i=i+1;
    dt2=ddt2*t+dt_2;
    t2=1/2*ddt2*t^2+dt_2*t+t_2;
syms t3 t4 t5
r__F =[l2*cos(t2) + l3*cos(t3) + l4*cos(t4);
 l2*sin(t2) + l3*sin(t3) + l4*sin(t4)];
r_F_real=[l1;0];
eqn1=r__F-r_F_real;
r__G =[l2*cos(t2) + l3*cos(t3) + l5*cos(t5) + lBE*cos(t4);
 l2*sin(t2) + l3*sin(t3) + l5*sin(t5) + lBE*sin(t4)];
r_G_2_real=h;
eqn2=-r__G(2)+r_G_2_real;
f=solve([eqn1;eqn2],[t3;t4;t5])
t3=f.t3(1);
t4=f.t4(1);
t5=f.t5(1);
if cos(t5)<0
t5=-(t5-pi);
end

q=mod([t3;t4;t5],2*pi);


syms dt3 dt4 dt5
v__F =[- dt2*l2*sin(t2) - dt3*l3*sin(t3) - dt4*l4*sin(t4);
   dt2*l2*cos(t2) + dt3*l3*cos(t3) + dt4*l4*cos(t4)]; %=0
   eqn1=v__F;
v__G =[- dt2*l2*sin(t2) - dt3*l3*sin(t3) - dt5*l5*sin(t5) - dt4*lBE*sin(t4);
   dt2*l2*cos(t2) + dt3*l3*cos(t3) + dt5*l5*cos(t5) + dt4*lBE*cos(t4)];
 eqn2=v__G(2);%=0
a=jacobian([eqn1;eqn2],[dt3;dt4;dt5]);
b=[eqn1;eqn2]-a*[dt3;dt4;dt5];
dq=-a\b;
dt3=dq(1);
dt4=dq(2);
dt5=dq(3);



MM =[ 0,        -l3*sin(t3),        -l4*sin(t4),                  0,              0,               0,              0,               0,              0,               0,               0,              0,                         0,                        0,              0,               0, 0
 0,         l3*cos(t3),         l4*cos(t4),                  0,              0,               0,              0,               0,              0,               0,               0,              0,                         0,                        0,              0,               0, 0
 0,         l3*cos(t3),        lBE*cos(t4),         l5*cos(t5),              0,               0,              0,               0,              0,               0,               0,              0,                         0,                        0,              0,               0, 0
 0,                  0,                  0,                  0,              1,               0,             -1,               0,              0,               0,               0,              0,                         0,                        0,              0,               0, 0
 0,                  0,                  0,                  0,              0,               1,              0,              -1,              0,               0,               0,              0,                         0,                        0,              0,               0, 0
 0,  (l3*m3*sin(t3))/2,                  0,                  0,              0,               0,              1,               0,             -1,               0,               0,              0,                         0,                        0,              0,               0, 0
 0, -(l3*m3*cos(t3))/2,                  0,                  0,              0,               0,              0,               1,              0,              -1,               0,              0,                         0,                        0,              0,               0, 0
 0,      l3*m4*sin(t3),  (l4*m4*sin(t4))/2,                  0,              0,               0,              0,               0,              1,               0,               1,              0,                        -1,                        0,              0,               0, 0
 0,     -l3*m4*cos(t3), -(l4*m4*cos(t4))/2,                  0,              0,               0,              0,               0,              0,               1,               0,              1,                         0,                       -1,              0,               0, 0
 0,      l3*m5*sin(t3),     lBE*m5*sin(t4),  (l5*m5*sin(t5))/2,              0,               0,              0,               0,              0,               0,               0,              0,                         1,                        0,             -1,               0, 0
 0,     -l3*m5*cos(t3),    -lBE*m5*cos(t4), -(l5*m5*cos(t5))/2,              0,               0,              0,               0,              0,               0,               0,              0,                         0,                        1,              0,              -1, 0
 0,      l3*m6*sin(t3),     lBE*m6*sin(t4),      l5*m6*sin(t5),              0,               0,              0,               0,              0,               0,               0,              0,                         0,                        0,              1,               0, 0
 0,     -l3*m6*cos(t3),    -lBE*m6*cos(t4),     -l5*m6*cos(t5),              0,               0,              0,               0,              0,               0,               0,              0,                         0,                        0,              0,               1, 1
 1,                  0,                  0,                  0, (l2*sin(t2))/2, -(l2*cos(t2))/2, (l2*sin(t2))/2, -(l2*cos(t2))/2,              0,               0,               0,              0,                         0,                        0,              0,               0, 0
 0,                -I3,                  0,                  0,              0,               0, (l3*sin(t3))/2, -(l3*cos(t3))/2, (l3*sin(t3))/2, -(l3*cos(t3))/2,               0,              0,                         0,                        0,              0,               0, 0
 0,                  0,                -I4,                  0,              0,               0,              0,               0, (l4*sin(t4))/2, -(l4*cos(t4))/2, -(l4*sin(t4))/2, (l4*cos(t4))/2, -(sin(t4)*(l4 - 2*lBE))/2, (cos(t4)*(l4 - 2*lBE))/2,              0,               0, 0
 0,                  0,                  0,                -I5,              0,               0,              0,               0,              0,               0,               0,              0,            (l5*sin(t5))/2,          -(l5*cos(t5))/2, (l5*sin(t5))/2, -(l5*cos(t5))/2, 0];
 
 
BB =[                                          - l2*cos(t2)*dt2^2 - l3*cos(t3)*dt3^2 - l4*cos(t4)*dt4^2 - ddt2*l2*sin(t2)
                                          - l2*sin(t2)*dt2^2 - l3*sin(t3)*dt3^2 - l4*sin(t4)*dt4^2 + ddt2*l2*cos(t2)
                      - l2*sin(t2)*dt2^2 - l3*sin(t3)*dt3^2 - lBE*sin(t4)*dt4^2 - l5*sin(t5)*dt5^2 + ddt2*l2*cos(t2)
                                                                            (l2*m2*(cos(t2)*dt2^2 + ddt2*sin(t2)))/2
                                                            - g*m2 - m2*((ddt2*l2*cos(t2))/2 - (dt2^2*l2*sin(t2))/2)
                                                  (m3*(2*l2*cos(t2)*dt2^2 + l3*cos(t3)*dt3^2 + 2*ddt2*l2*sin(t2)))/2
                                         -(m3*(- 2*l2*sin(t2)*dt2^2 - l3*sin(t3)*dt3^2 + 2*g + 2*ddt2*l2*cos(t2)))/2
                             (m4*(2*l2*cos(t2)*dt2^2 + 2*l3*cos(t3)*dt3^2 + l4*cos(t4)*dt4^2 + 2*ddt2*l2*sin(t2)))/2
                       (m4*(2*l2*sin(t2)*dt2^2 + 2*l3*sin(t3)*dt3^2 + l4*sin(t4)*dt4^2 - 2*g - 2*ddt2*l2*cos(t2)))/2
       (m5*(2*l2*cos(t2)*dt2^2 + 2*l3*cos(t3)*dt3^2 + 2*lBE*cos(t4)*dt4^2 + l5*cos(t5)*dt5^2 + 2*ddt2*l2*sin(t2)))/2
 (m5*(2*l2*sin(t2)*dt2^2 + 2*l3*sin(t3)*dt3^2 + 2*lBE*sin(t4)*dt4^2 + l5*sin(t5)*dt5^2 - 2*g - 2*ddt2*l2*cos(t2)))/2
                   m6*(l2*cos(t2)*dt2^2 + l3*cos(t3)*dt3^2 + lBE*cos(t4)*dt4^2 + l5*cos(t5)*dt5^2 + ddt2*l2*sin(t2))
               m6*(l2*sin(t2)*dt2^2 + l3*sin(t3)*dt3^2 + lBE*sin(t4)*dt4^2 + l5*sin(t5)*dt5^2 - g - ddt2*l2*cos(t2))
                                                                                                            -I2*ddt2
                                                                                                                   0
                                                                                                                   0
                                                                                                                   0];
un=-MM\BB;
     T(i)=un(1);
     %For calculation
     ddt3=un(2);
     ddt4=un(3);
     ddt5=un(4);
 F_D_x(i)=un(5);
 F_D_y(i)=un(6);
 F_A_x(i)=un(7);
 F_A_y(i)=un(8);
 F_B_x(i)=un(9);
 F_B_y(i)=un(10);
 F_F_x(i)=un(11);
 F_F_y(i)=un(12);
 F_E_x(i)=un(13);
 F_E_y(i)=un(14);
 F_G_x(i)=un(15);
 F_G_y(i)=un(16);
  F_n(i)=un(17);
  %For ploting
 ddt_2(i)=ddt2;
     ddt_3(i)=un(2);
     ddt_4(i)=un(3);
     ddt_5(i)=un(4);
  
  r_A(:,i) =[l2*cos(t2)
 l2*sin(t2)];
 
r_B(:,i) =[l2*cos(t2) + l3*cos(t3)
 l2*sin(t2) + l3*sin(t3)];
 
 
r_F(:,i) =[l2*cos(t2) + l3*cos(t3) + l4*cos(t4)
 l2*sin(t2) + l3*sin(t3) + l4*sin(t4)];
 
 
r_E(:,i) =[l2*cos(t2) + l3*cos(t3) + lBE*cos(t4)
 l2*sin(t2) + l3*sin(t3) + lBE*sin(t4)];
 
 
r_G(:,i) =[l2*cos(t2) + l3*cos(t3) + l5*cos(t5) + lBE*cos(t4)
 l2*sin(t2) + l3*sin(t3) + l5*sin(t5) + lBE*sin(t4)];
 
 
r_X(:,i) =[X*cos(t5) + l2*cos(t2) + l3*cos(t3) + lBE*cos(t4)
 X*sin(t5) + l2*sin(t2) + l3*sin(t3) + lBE*sin(t4)];
 
 
r_g_2(:,i) =[(l2*cos(t2))/2
 (l2*sin(t2))/2];
 
 
r_g_3(:,i) =[l2*cos(t2) + (l3*cos(t3))/2
 l2*sin(t2) + (l3*sin(t3))/2];
 
 
r_g_4(:,i) =[l2*cos(t2) + l3*cos(t3) + (l4*cos(t4))/2
 l2*sin(t2) + l3*sin(t3) + (l4*sin(t4))/2];
 
 
r_g_5(:,i) =[l2*cos(t2) + l3*cos(t3) + (l5*cos(t5))/2 + lBE*cos(t4)
 l2*sin(t2) + l3*sin(t3) + (l5*sin(t5))/2 + lBE*sin(t4)];
 
 
v_A(:,i) =[-dt2*l2*sin(t2)
  dt2*l2*cos(t2)];
 
 
v_B(:,i) =[- dt2*l2*sin(t2) - dt3*l3*sin(t3)
   dt2*l2*cos(t2) + dt3*l3*cos(t3)];
 
 
v_F(:,i) =[- dt2*l2*sin(t2) - dt3*l3*sin(t3) - dt4*l4*sin(t4)
   dt2*l2*cos(t2) + dt3*l3*cos(t3) + dt4*l4*cos(t4)];
 
 
v_E(:,i) =[- dt2*l2*sin(t2) - dt3*l3*sin(t3) - dt4*lBE*sin(t4)
   dt2*l2*cos(t2) + dt3*l3*cos(t3) + dt4*lBE*cos(t4)];
 
 
v_G(:,i) =[- dt2*l2*sin(t2) - dt3*l3*sin(t3) - dt5*l5*sin(t5) - dt4*lBE*sin(t4)
   dt2*l2*cos(t2) + dt3*l3*cos(t3) + dt5*l5*cos(t5) + dt4*lBE*cos(t4)];
 
 
v_X(:,i) =[- X*dt5*sin(t5) - dt2*l2*sin(t2) - dt3*l3*sin(t3) - dt4*lBE*sin(t4)
   dt2*l2*cos(t2) + dt3*l3*cos(t3) + dt4*lBE*cos(t4) + X*dt5*cos(t5)];
 
 
v_g_2(:,i) =[-(dt2*l2*sin(t2))/2
  (dt2*l2*cos(t2))/2];
 
 
v_g_3(:,i) =[- dt2*l2*sin(t2) - (dt3*l3*sin(t3))/2
   dt2*l2*cos(t2) + (dt3*l3*cos(t3))/2];
 
 
v_g_4(:,i) =[- dt2*l2*sin(t2) - dt3*l3*sin(t3) - (dt4*l4*sin(t4))/2
   dt2*l2*cos(t2) + dt3*l3*cos(t3) + (dt4*l4*cos(t4))/2];
 
 
v_g_5(:,i) =[- dt2*l2*sin(t2) - dt3*l3*sin(t3) - (dt5*l5*sin(t5))/2 - dt4*lBE*sin(t4)
   dt2*l2*cos(t2) + dt3*l3*cos(t3) + (dt5*l5*cos(t5))/2 + dt4*lBE*cos(t4)];
 
 
a_A(:,i) =[- l2*cos(t2)*dt2^2 - ddt2*l2*sin(t2)
 - l2*sin(t2)*dt2^2 + ddt2*l2*cos(t2)];
 
 
a_B(:,i) =[- l2*cos(t2)*dt2^2 - l3*cos(t3)*dt3^2 - ddt2*l2*sin(t2) - ddt3*l3*sin(t3)
 - l2*sin(t2)*dt2^2 - l3*sin(t3)*dt3^2 + ddt2*l2*cos(t2) + ddt3*l3*cos(t3)];
 
 
a_F(:,i) =[- l2*cos(t2)*dt2^2 - l3*cos(t3)*dt3^2 - l4*cos(t4)*dt4^2 - ddt2*l2*sin(t2) - ddt3*l3*sin(t3) - ddt4*l4*sin(t4)
 - l2*sin(t2)*dt2^2 - l3*sin(t3)*dt3^2 - l4*sin(t4)*dt4^2 + ddt2*l2*cos(t2) + ddt3*l3*cos(t3) + ddt4*l4*cos(t4)];
 
 
a_E(:,i) =[- l2*cos(t2)*dt2^2 - l3*cos(t3)*dt3^2 - lBE*cos(t4)*dt4^2 - ddt2*l2*sin(t2) - ddt3*l3*sin(t3) - ddt4*lBE*sin(t4)
 - l2*sin(t2)*dt2^2 - l3*sin(t3)*dt3^2 - lBE*sin(t4)*dt4^2 + ddt2*l2*cos(t2) + ddt3*l3*cos(t3) + ddt4*lBE*cos(t4)];
 
 
a_G(:,i) =[- l2*cos(t2)*dt2^2 - l3*cos(t3)*dt3^2 - lBE*cos(t4)*dt4^2 - l5*cos(t5)*dt5^2 - ddt2*l2*sin(t2) - ddt3*l3*sin(t3) - ddt5*l5*sin(t5) - ddt4*lBE*sin(t4)
 - l2*sin(t2)*dt2^2 - l3*sin(t3)*dt3^2 - lBE*sin(t4)*dt4^2 - l5*sin(t5)*dt5^2 + ddt2*l2*cos(t2) + ddt3*l3*cos(t3) + ddt5*l5*cos(t5) + ddt4*lBE*cos(t4)];
 
 
a_X(:,i) =[- l2*cos(t2)*dt2^2 - l3*cos(t3)*dt3^2 - lBE*cos(t4)*dt4^2 - X*cos(t5)*dt5^2 - X*ddt5*sin(t5) - ddt2*l2*sin(t2) - ddt3*l3*sin(t3) - ddt4*lBE*sin(t4)
 - l2*sin(t2)*dt2^2 - l3*sin(t3)*dt3^2 - lBE*sin(t4)*dt4^2 - X*sin(t5)*dt5^2 + ddt2*l2*cos(t2) + ddt3*l3*cos(t3) + ddt4*lBE*cos(t4) + X*ddt5*cos(t5)];
 
 
a_g_2(:,i) =[- (ddt2*l2*sin(t2))/2 - (dt2^2*l2*cos(t2))/2
   (ddt2*l2*cos(t2))/2 - (dt2^2*l2*sin(t2))/2];
 
 
a_g_3(:,i) =[- l2*cos(t2)*dt2^2 - ddt2*l2*sin(t2) - (ddt3*l3*sin(t3))/2 - (dt3^2*l3*cos(t3))/2
 - l2*sin(t2)*dt2^2 + ddt2*l2*cos(t2) + (ddt3*l3*cos(t3))/2 - (dt3^2*l3*sin(t3))/2];
 
 
a_g_4(:,i) =[- l2*cos(t2)*dt2^2 - l3*cos(t3)*dt3^2 - ddt2*l2*sin(t2) - ddt3*l3*sin(t3) - (ddt4*l4*sin(t4))/2 - (dt4^2*l4*cos(t4))/2
 - l2*sin(t2)*dt2^2 - l3*sin(t3)*dt3^2 + ddt2*l2*cos(t2) + ddt3*l3*cos(t3) + (ddt4*l4*cos(t4))/2 - (dt4^2*l4*sin(t4))/2];
 
 
a_g_5(:,i) =[- l2*cos(t2)*dt2^2 - l3*cos(t3)*dt3^2 - lBE*cos(t4)*dt4^2 - ddt2*l2*sin(t2) - ddt3*l3*sin(t3) - (ddt5*l5*sin(t5))/2 - ddt4*lBE*sin(t4) - (dt5^2*l5*cos(t5))/2
 - l2*sin(t2)*dt2^2 - l3*sin(t3)*dt3^2 - lBE*sin(t4)*dt4^2 + ddt2*l2*cos(t2) + ddt3*l3*cos(t3) + (ddt5*l5*cos(t5))/2 + ddt4*lBE*cos(t4) - (dt5^2*l5*sin(t5))/2];
 
z(i,:)=[t2,q.',dt2,dq.',ddt2,ddt3,ddt4,ddt5];
end %calculate for plot

t=[0:dt:time];


% animation_1(t,z)

for i=1:1 %plot
figure  %theta
subplot(2,2,1)
plot(t,z(:,1))
xlabel('Time (S)')
ylabel('angle of link 2 (rad)')
subplot(2,2,2)
plot(t,z(:,2))
xlabel('Time (S)')
ylabel('angle of link 3 (rad)')
subplot(2,2,3)
plot(t,z(:,3))
xlabel('Time (S)')
ylabel('angle of link 4 (rad)')
subplot(2,2,4)
plot(t,z(:,4))
xlabel('Time (S)')
ylabel('angle of link 5 (rad)')
suptitle('Angle')


figure  %omega
subplot(2,2,1)
plot(t,z(:,5))
xlabel('Time (S)')
ylabel('omega of link 2 (rad/s)')
subplot(2,2,2)
plot(t,z(:,6))
xlabel('Time (S)')
ylabel('omega of link 3 (rad/s)')
subplot(2,2,3)
plot(t,z(:,7))
xlabel('Time (S)')
ylabel('omega of link 4 (rad/s)')
subplot(2,2,4)
plot(t,z(:,8))
xlabel('Time (S)')
ylabel('omega of link 5 (rad/s)')
suptitle('Omega')

figure  %alpha
subplot(2,2,1)
plot(t,ddt_2)
xlabel('Time (S)')
ylabel('alpha of link 2 (rad/s^2)')
subplot(2,2,2)
plot(t,ddt_3)
xlabel('Time (S)')
ylabel('alpha of link 3 (rad/s^2)')
subplot(2,2,3)
plot(t,ddt_4)
xlabel('Time (S)')
ylabel('alpha of link 4 (rad/s^2)')
subplot(2,2,4)
plot(t,ddt_5)
xlabel('Time (S)')
ylabel('alpha of link 5 (rad/s^2)')
suptitle('Alpha')

figure %joint D
subplot(1,2,1)
plot(t,F_D_x)
xlabel('Time (S)')
ylabel('F_D_x (N)')
subplot(1,2,2)
plot(t,F_D_y)
xlabel('Time (S)')
ylabel('F_D_y (N)')
suptitle('Joint D')

figure %joint F
subplot(1,2,1)
plot(t,F_F_x)
xlabel('Time (S)')
ylabel('F_F_x (N)')
subplot(1,2,2)
plot(t,F_F_y)
xlabel('Time (S)')
ylabel('F_F_y (N)')
suptitle('Joint F')

figure %joint N
plot(t,F_n)
xlabel('Time (S)')
ylabel('F_N_Y (N)')
suptitle('Joint N')


figure %joint A
subplot(4,2,1)
plot(t,r_A(1,:))
xlabel('Time (S)')
ylabel('X_A (m)')
subplot(4,2,2)
plot(t,r_A(2,:))
xlabel('Time (S)')
ylabel('Y_A (m)')
subplot(4,2,3)
plot(t,v_A(1,:))
xlabel('Time (S)')
ylabel('v_A_X (m/s)')
subplot(4,2,4)
plot(t,v_A(2,:))
xlabel('Time (S)')
ylabel('v_A_Y (m/s)')
subplot(4,2,5)
plot(t,a_A(1,:))
xlabel('Time (S)')
ylabel('a_A_X (m/s^2)')
subplot(4,2,6)
plot(t,a_A(2,:))
xlabel('Time (S)')
ylabel('a_A_Y (m/s^2)')
subplot(4,2,7)


subplot(4,2,8)
plot(t,F_A_y)
xlabel('Time (S)')
ylabel('F_A_Y (N)')
suptitle('Joint A')

figure %joint B
subplot(4,2,1)
plot(t,r_B(1,:))
xlabel('Time (S)')
ylabel('X_B (m)')
subplot(4,2,2)
plot(t,r_B(2,:))
xlabel('Time (S)')
ylabel('Y_B (m)')
subplot(4,2,3)
plot(t,v_B(1,:))
xlabel('Time (S)')
ylabel('v_B_X (m/s)')
subplot(4,2,4)
plot(t,v_B(2,:))
xlabel('Time (S)')
ylabel('v_B_Y (m/s)')
subplot(4,2,5)
plot(t,a_B(1,:))
xlabel('Time (S)')
ylabel('a_B_X (m/s^2)')
subplot(4,2,6)
plot(t,a_B(2,:))
xlabel('Time (S)')
ylabel('a_B_Y (m/s^2)')
subplot(4,2,7)
plot(t,F_B_x)
xlabel('Time (S)')
ylabel('F_B_X (N)')
subplot(4,2,8)
plot(t,F_B_y)
xlabel('Time (S)')
ylabel('F_B_Y (N)')
suptitle('Joint B')

figure %joint E
subplot(4,2,1)
plot(t,r_E(1,:))
xlabel('Time (S)')
ylabel('X_E (m)')
subplot(4,2,2)
plot(t,r_E(2,:))
xlabel('Time (S)')
ylabel('Y_E (m)')
subplot(4,2,3)
plot(t,v_E(1,:))
xlabel('Time (S)')
ylabel('v_E_X (m/s)')
subplot(4,2,4)
plot(t,v_E(2,:))
xlabel('Time (S)')
ylabel('v_E_Y (m/s)')
subplot(4,2,5)
plot(t,a_E(1,:))
xlabel('Time (S)')
ylabel('a_E_X (m/s^2)')
subplot(4,2,6)
plot(t,a_E(2,:))
xlabel('Time (S)')
ylabel('a_E_Y (m/s^2)')
subplot(4,2,7)
plot(t,F_E_x)
xlabel('Time (S)')
ylabel('F_E_X (N)')
subplot(4,2,8)
plot(t,F_E_y)
xlabel('Time (S)')
ylabel('F_E_Y (N)')
suptitle('Joint E')

figure %joint G
subplot(4,2,1)
plot(t,r_G(1,:))
xlabel('Time (S)')
ylabel('X_G (m)')
subplot(4,2,2)
plot(t,r_G(2,:))
xlabel('Time (S)')
ylabel('Y_G (m)')
subplot(4,2,3)
plot(t,v_G(1,:))
xlabel('Time (S)')
ylabel('v_G_X (m/s)')
subplot(4,2,4)
plot(t,v_G(2,:))
xlabel('Time (S)')
ylabel('v_G_Y (m/s)')
subplot(4,2,5)
plot(t,a_G(1,:))
xlabel('Time (S)')
ylabel('a_G_X (m/s^2)')
subplot(4,2,6)
plot(t,a_G(2,:))
xlabel('Time (S)')
ylabel('a_G_Y (m/s^2)')
subplot(4,2,7)
plot(t,F_G_x)
xlabel('Time (S)')
ylabel('F_G_X (N)')
subplot(4,2,8)
plot(t,F_G_y)
xlabel('Time (S)')
ylabel('F_G_Y (N)')
suptitle('Joint G')

figure %point C_2
subplot(3,2,1)
plot(t,r_g_2(1,:))
xlabel('Time (S)')
ylabel('X_C_2 (m)')
subplot(3,2,2)
plot(t,r_g_2(2,:))
xlabel('Time (S)')
ylabel('Y_C_2 (m)')
subplot(3,2,3)
plot(t,v_g_2(1,:))
xlabel('Time (S)')
ylabel('v_C_x_2 (m/s)')
subplot(3,2,4)
plot(t,v_g_2(2,:))
xlabel('Time (S)')
ylabel('v_C_Y_2 (m/s)')
subplot(3,2,5)
plot(t,a_g_2(1,:))
xlabel('Time (S)')
ylabel('a_C_X_2 (m/s^2)')
subplot(3,2,6)
plot(t,a_g_2(2,:))
xlabel('Time (S)')
ylabel('a_C_Y_2 (m/s^2)')
suptitle('C_2')

figure %point C_3
subplot(3,2,1)
plot(t,r_g_3(1,:))
xlabel('Time (S)')
ylabel('X_C_3 (m)')
subplot(3,2,2)
plot(t,r_g_3(2,:))
xlabel('Time (S)')
ylabel('Y_C_3 (m)')
subplot(3,2,3)
plot(t,v_g_3(1,:))
xlabel('Time (S)')
ylabel('v_C_x_3 (m/s)')
subplot(3,2,4)
plot(t,v_g_3(2,:))
xlabel('Time (S)')
ylabel('v_C_Y_3 (m/s)')
subplot(3,2,5)
plot(t,a_g_3(1,:))
xlabel('Time (S)')
ylabel('a_C_X_3 (m/s^2)')
subplot(3,2,6)
plot(t,a_g_3(2,:))
xlabel('Time (S)')
ylabel('a_C_Y_3 (m/s^2)')
suptitle('C_3')

figure %point C_4
subplot(3,2,1)
plot(t,r_g_4(1,:))
xlabel('Time (S)')
ylabel('X_C_4 (m)')
subplot(3,2,2)
plot(t,r_g_4(2,:))
xlabel('Time (S)')
ylabel('Y_C_4 (m)')
subplot(3,2,3)
plot(t,v_g_4(1,:))
xlabel('Time (S)')
ylabel('v_C_x_4 (m/s)')
subplot(3,2,4)
plot(t,v_g_4(2,:))
xlabel('Time (S)')
ylabel('v_C_Y_4 (m/s)')
subplot(3,2,5)
plot(t,a_g_4(1,:))
xlabel('Time (S)')
ylabel('a_C_X_4 (m/s^2)')
subplot(3,2,6)
plot(t,a_g_4(2,:))
xlabel('Time (S)')
ylabel('a_C_Y_4 (m/s^2)')
suptitle('C_4')

figure %point C_5
subplot(3,2,1)
plot(t,r_g_5(1,:))
xlabel('Time (S)')
ylabel('X_C_5 (m)')
subplot(3,2,2)
plot(t,r_g_5(2,:))
xlabel('Time (S)')
ylabel('Y_C_5 (m)')
subplot(3,2,3)
plot(t,v_g_5(1,:))
xlabel('Time (S)')
ylabel('v_C_x_5 (m/s)')
subplot(3,2,4)
plot(t,v_g_5(2,:))
xlabel('Time (S)')
ylabel('v_C_Y_5 (m/s)')
subplot(3,2,5)
plot(t,a_g_5(1,:))
xlabel('Time (S)')
ylabel('a_C_X_5 (m/s^2)')
subplot(3,2,6)
plot(t,a_g_5(2,:))
xlabel('Time (S)')
ylabel('a_C_Y_5 (m/s^2)')
suptitle('C_5')

figure %point X
subplot(3,2,1)
plot(t,r_X(1,:))
xlabel('Time (S)')
ylabel('X (m)')
subplot(3,2,2)
plot(t,r_X(2,:))
xlabel('Time (S)')
ylabel('Y (m)')
subplot(3,2,3)
plot(t,v_X(1,:))
xlabel('Time (S)')
ylabel('v_X (m/s)')
subplot(3,2,4)
plot(t,v_X(2,:))
xlabel('Time (S)')
ylabel('v_Y (m/s)')
subplot(3,2,5)
plot(t,a_X(1,:))
xlabel('Time (S)')
ylabel('a_X (m/s^2)')
subplot(3,2,6)
plot(t,a_X(2,:))
xlabel('Time (S)')
ylabel('a_Y (m/s^2)')
suptitle('Point X on EG')

figure %tourqe T
plot(t,T)
xlabel('Time (S)')
ylabel('Tourqe link 2')
suptitle('Tourqe T')
end

%{
function a=animation(t,z)
a=0;
global l1 l2 l3 l4 l5 lBE h X T m2 m3 m4 m5 m6 I2 I3 I4 I5 g
for i=1:length(t)
      t2=z(i,1);
      t3=z(i,2);
      t4=z(i,3);
      t5=z(i,4);
     
  r_A =[l2*cos(t2)
 l2*sin(t2)];
 
r_B =[l2*cos(t2) + l3*cos(t3)
 l2*sin(t2) + l3*sin(t3)];
 
r_F =[l2*cos(t2) + l3*cos(t3) + l4*cos(t4)
 l2*sin(t2) + l3*sin(t3) + l4*sin(t4)];
 
r_E =[l2*cos(t2) + l3*cos(t3) + lBE*cos(t4)
 l2*sin(t2) + l3*sin(t3) + lBE*sin(t4)];
 
r_G =[l2*cos(t2) + l3*cos(t3) + l5*cos(t5) + lBE*cos(t4)
 l2*sin(t2) + l3*sin(t3) + l5*sin(t5) + lBE*sin(t4)];
 
 plot([0 r_A(1)],[0 r_A(2)])
 hold on 
 plot([r_B(1) r_A(1)],[r_B(2) r_A(2)])
 plot([r_B(1) r_F(1)],[r_B(2) r_F(2)])
 plot([r_E(1) r_G(1)],[r_E(2) r_G(2)])
 xlim([-3 14])
 ylim([-3 8])
 theta=0:.1:2.1*pi;
 plot(.5*cos(theta),.5*sin(theta))
 plot(l1+.5*cos(theta),.5*sin(theta))
 plot(r_A(1)+.2*cos(theta),r_A(2)+.2*sin(theta))
 plot(r_B(1)+.2*cos(theta),r_B(2)+.2*sin(theta))
 w=.4;
 b=.8;
 plot(r_G(1)+[-b b b -b -b],r_G(2)+[-w -w w w -w])
 xlabel('X (m)')
ylabel('Y (m)')
title(['Time= ',num2str(i*(t(2)-t(1))),'s'])
 hold off
 pause(.1)
end
end 
%}