hold off 
t=0;
% for t=0:0.1:1;
    x1=0:0.01:1 ;
    Nt =  1 + floor(10*abs(sin(0.5*pi*t)));
    max  = (1/(2.0*Nt)+0.1)*sin(2*Nt*pi*x1);
    f1 = x1;
    f2 = 1- f1;
    plot(f1,f2,'Linewidth',2)
    hold on
    pause(0.1);
 %end;
hold on
A=importdata('.\pf_DF9_0.dat');
F1=A(:,1)
F2=A(:,2)
plot(F1,F2,'ro'); %»­¸öÌå
hold on
