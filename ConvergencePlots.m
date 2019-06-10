
%type = 'No Interface';
type = 'One Interface';
switch type

    case 'No Interface'

h     = [0.1,0.05,0.025,0.0125,0.00625,0.003125,...
         0.0015625, 0.00078125];

Elerr = [0.3235369176998824,0.22304962790032315,...
         0.13181185314682203,0.07132186423285303,...
         0.03700832358191203,0.018836424316228784,...
         0.009500465043459295,0.004770679242519312];
         


    
Hferr = [0.3303304317175603,0.224501282835075,...
         0.13234322102634932,0.07147950975215762,...
         0.03705058312518213,0.01884731154135237,...
         0.009503224067225651,0.004771373368334899];
     
     
     

figure(1)

clf

loglog(h,Elerr,'o')
%plot(log(h),log(ElectricError));

hold on


loglog(h,Elerr)
%loglog(h,10^(1)*h.^(1.8))
hold on
%Pick a basis point for the triangle
xseed = 0.001;
yseed = 0.007;
%desiredSlope Of triangle
slope = 1;
%Another x point
%xnext = h(4);
xnext = h(6);

b = log10(yseed/(xseed^slope));


ynext = 10^(b)*xnext^slope;


x = [xseed, xseed, xnext, xseed];
y = [yseed, ynext, ynext, yseed];

loglog(x,y,'r')

text(xseed*0.7+0.3*xnext,ynext+0.001+0.001,'1')
text(xseed-0.0001,0.5*yseed+0.5*ynext,'1')





title('L2 Error on Electric Field Vs. Mesh Size')

%Now the ticks
%The x-axis
leftaxis = 0.0006;
rightaxis = 0.1+0.01;
upaxis = 0.4;
downaxis = 0.0043;

axis([leftaxis rightaxis downaxis upaxis])
xlabel('Mesh Size')

ax = log10(leftaxis);
bx = log10(rightaxis);
tx = ax:(bx-ax)/10:bx;
tx = 10.^(tx);
ay = log10(downaxis);
by = log10(upaxis);
ty = ay:(by-ay)/10:by;
ty = 10.^(ty);

ylabel('L2 Error')

xticks(tx)
yticks(ty)
grid on

 figure(2)

clf

loglog(h,Hferr,'o')
%plot(log(h),log(ElectricError));

hold on


loglog(h,Hferr)
%loglog(h,10^(1)*h.^(1.8))
hold on
%Pick a basis point for the triangle
xseed = 0.001;
yseed = 0.007;
%desiredSlope Of triangle
slope = 1;
%Another x point
%xnext = h(4);
xnext = h(6);

b = log10(yseed/(xseed^slope));


ynext = 10^(b)*xnext^slope;


x = [xseed, xseed, xnext, xseed];
y = [yseed, ynext, ynext, yseed];

loglog(x,y,'r')

text(xseed*0.7+0.3*xnext,ynext+0.001+0.001,'1')
text(xseed-0.0001,0.5*yseed+0.5*ynext,'1')





title('L2 Error on Magnetic Field Vs. Mesh Size')

%Now the ticks
%The x-axis
leftaxis = 0.0006;
rightaxis = 0.1+0.01;
upaxis = 0.4;
downaxis = 0.0043;

axis([leftaxis rightaxis downaxis upaxis])
xlabel('Mesh Size')

ax = log10(leftaxis);
bx = log10(rightaxis);
tx = ax:(bx-ax)/10:bx;
tx = 10.^(tx);
ay = log10(downaxis);
by = log10(upaxis);
ty = ay:(by-ay)/10:by;
ty = 10.^(ty);

ylabel('L2 Error')

xticks(tx)
yticks(ty)
grid on
     
    case 'One Interface'

   h  = [0.1,0.05,0.025,0.0125,0.00625,0.003125,...
         0.0015625, 0.00078125];

Elerr = [1.5020192828751187,0.32055092292529375,...
         0.04804850245582069, 0.010463434940433185,...
         0.0024432682607272507,0.0005855334883575576,...
         0.00014367790574039885,3.557675098710133e-05];
         


    
Hferr = [1.062725826917307,0.3558644362937835,...
         0.07188692174401906,0.0156968014580204,...
         0.00369170129490299,0.0008962230813669055,...
         0.00022071520959128732,5.476260090127839e-05];
     
     
     

figure(1)

clf

loglog(h,Elerr,'o')
%plot(log(h),log(ElectricError));

hold on


loglog(h,Elerr)
%loglog(h,10^(1)*h.^(1.8))
hold on
%Pick a basis point for the triangle
xseed = 0.001;
yseed = 0.0001;
%desiredSlope Of triangle
slope = 2;
%Another x point
%xnext = h(4);
xnext = h(6);

b = log10(yseed/(xseed^slope));


ynext = 10^(b)*xnext^slope;


x = [xseed, xseed, xnext, xseed];
y = [yseed, ynext, ynext, yseed];

loglog(x,y,'r')

text(xseed*0.7+0.3*xnext,ynext+0.00025,'1')
text(xseed-0.0001,0.5*yseed+0.5*ynext-0.0002,'2')





title('L2 Error on Electric Field Vs. Mesh Size')

%Now the ticks
%The x-axis
leftaxis = 0.0006;
rightaxis = 0.1+0.01;
upaxis = 1.8;
downaxis = 0.00003;

axis([leftaxis rightaxis downaxis upaxis])
xlabel('Mesh Size')

ax = log10(leftaxis);
bx = log10(rightaxis);
tx = ax:(bx-ax)/10:bx;
tx = 10.^(tx);
ay = log10(downaxis);
by = log10(upaxis);
ty = ay:(by-ay)/10:by;
ty = 10.^(ty);

ylabel('L2 Error')

xticks(tx)
yticks(ty)
grid on

 figure(2)

clf

loglog(h,Hferr,'o')
%plot(log(h),log(ElectricError));

hold on


loglog(h,Hferr)
%loglog(h,10^(1)*h.^(1.8))
hold on
%Pick a basis point for the triangle
xseed = 0.001;
yseed = 0.0002;
%desiredSlope Of triangle
slope = 2;
%Another x point
%xnext = h(4);
xnext = h(6);

b = log10(yseed/(xseed^slope));


ynext = 10^(b)*xnext^slope;


x = [xseed, xseed, xnext, xseed];
y = [yseed, ynext, ynext, yseed];

loglog(x,y,'r')

text(xseed*0.7+0.3*xnext,ynext+0.0005,'1')
text(xseed-0.0001,0.5*yseed+0.5*ynext-0.0003,'2')





title('L2 Error on Magnetic Field Vs. Mesh Size')

%Now the ticks
%The x-axis
leftaxis = 0.0006;
rightaxis = 0.1+0.01;
upaxis = 1.8;
downaxis = 0.00003;

axis([leftaxis rightaxis downaxis upaxis])
xlabel('Mesh Size')

ax = log10(leftaxis);
bx = log10(rightaxis);
tx = ax:(bx-ax)/10:bx;
tx = 10.^(tx);
ay = log10(downaxis);
by = log10(upaxis);
ty = ay:(by-ay)/10:by;
ty = 10.^(ty);

ylabel('L2 Error')

xticks(tx)
yticks(ty)
grid on
    
end   
 
