function a = animation_1(t,z)
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
 xlim([-l2 l1+l2+.05])
 ylim([-l2 h+l2/2])
 theta=0:.1:2.1*pi;
 %circle
 plot(1*cos(theta),1*sin(theta))
 plot(l1+1*cos(theta),1*sin(theta))
 plot(r_A(1)+.008*cos(theta),r_A(2)+.008*sin(theta))
 plot(r_B(1)+.008*cos(theta),r_B(2)+.008*sin(theta))
 plot(r_E(1)+.008*cos(theta),r_E(2)+.008*sin(theta))
 %rectangle 
 w=.008;
 b=.016;
 plot(r_G(1)+[-b b b -b -b],r_G(2)+[-w -w w w -w])
 xlabel('X (m)')
ylabel('Y (m)')
title(['Time= ',num2str(i*(t(2)-t(1))),'s'])
 hold off
 pause(.05)
end
end