m=15.2564;
g=9.81;
w=m*g;
V=2.378;
B=.254;
Clb=w/(.5*998*V^2*B^2);
Cl0=0;
test=.9446-.0065*19.5*.9446^.6
       
Cl0=0;
i=1;

while i>1E-6;
    Clb=w/(.5*998*V^2*B^2);
    Cl0=Cl0+.0001;
    a=Cl0;
    b=a-.0065*19.5*a^.6;
    i=Clb-b;
end


M=0;
a=.0778-.25*B*tan(19.5*pi/180);
nu=8.721*10^(-7);
DCf=0.0004;
f=0;
z=1;
Fb=1.5
t=8;
while abs(z) > 1E-1
    t=t+.01;
    lambdaw=2;
    oldb=0;
    for i=2:.001:4
        b=t^1.1*(0.0120*(i^(.5))+0.0055*i^(5/2)/(Fb^2))
        if abs(b-Cl0)< abs(oldb-Cl0)
            lambdaw=i;
        end
        oldb=b;
    end 
t=t*pi/180;
Lm=lambdaw*B;
Lk=Lm+B*tan(19.5*pi/180)/(2*pi*tan(t));
Lc=2*Lm-Lk;
Cld=0.0120*lambdaw^.5*(t*180/pi)^1.1-0.0065*19.5*(0.0120*lambdaw^.5*(t*180/pi)^1.1)^.6;
Vm=V*(1-Cld/(lambdaw*cos(t)))^.5;
h=Lk*sin(t);
Lp=(.75-1/(5.21*(Fb^2/lambdaw^2)+2.39))*Lm;
c=.391-Lp;
Rn=(Vm*Lm)/nu;
Cf=0.075/((log10(Rn)-2)^2);
Sf=(lambdaw*B^2)/cos(19.5*pi/180);
Rv=.5*998*Vm^2*Sf*(Cf+DCf);
Rt=w*tan(t)+Rv/cos(t);
Mt=w*(c/cos(t)-f*sin(t))+Rv*(a-f);
z=M-Mt;
t=180/pi*t;
end