clc;
clear all;
close all;
x=input('Enter the input sequence');
N=length(x);
m=0:N-1;
X=zeros(1,N);
for k=1:N
    for n=1:N
        X(k)=X(k)+x(n)exp(-i*2*pi(k-1)*(n-1)/N);
    end
end

subplot(3,2,1);
stem(m,real(x));
xlabel('n');
ylabel('real x[n]');
title('REAL PART OF x[n]');
subplot(3,2,2);
stem(m,imag(x));
xlabel('n');
ylabel('imaginary x[n]');
title('IMAGINARY PART OF x[n]');
disp(X)


y=zeros(1,N);
for n=1:N
    for k=1:N
        y(n)=y(n)+X(k)exp(i*2*pi(k-1)*(n-1)/N);
    end
    y(n)=y(n)/N;
end
subplot(3,2,3);
stem(m,real(X));
xlabel('k');
ylabel('realX[k]');
title('REAL PART OF DFT');
subplot(3,2,4);
stem(m,imag(X));
xlabel('k');
ylabel('imaginaryX[k]');
title('IMAGINARY PART OF DFT')
subplot(3,2,5);
stem(m,real(y));
xlabel('n');
ylabel('real y[n]');
title(' REAL PART OF IDFT of X');
subplot(3,2,6);
stem(m,imag(y));
xlabel('n');
ylabel('imaginary y[n]');
title(' IMAGINARY PART OF IDFT of X')
disp(y)

x=[1,-1,1,-1,1,2,1,1];
l=length(x);
m=l/2;
n=l/4;
x1=zeros(1,m);
x2=zeros(1,m);
for a=1:4
    b=exp(-1i*2*pi*(a-1)/8);
    x1(a)=x(a)+x(a+4);
    x2(a)=b.*(x(a)-x(a+4));
    
end
disp(x1);
disp(x2);
x3=zeros(1,n);
x4=zeros(1,n);
x5=zeros(1,n);
x6=zeros(1,n);
for c=1:2
    x3(c)=x1(c)+x1(c+2);
    x4(c)=(-1i).*(x1(c)-x1(c+2));
    x5(c)=x2(c)+x2(c+2);
    x6(c)=(-1i).*(x2(c)-x2(c+2));
end
disp(x3);
disp(x4);
disp(x5);
disp(x6);
y=zeros(1,l);
y(1)=x3(1)+x3(2);
y(5)=x3(1)-x3(2);
y(2)=x5(1)+x5(2);
y(6)=x5(1)-x5(2);
y(3)=x4(1)+x4(2);
y(7)=x4(1)-x4(2);
y(4)=x6(1)+x6(2);
y(8)=x6(1)-x6(2);
disp(y)

DIT
x=[1,2,3,4,4,3,2,1];
l=length(x);
m=l/4;
n=l/2;
x1=bitrevorder(x);
x2=zeros(1,m);
x3=zeros(1,m);
x4=zeros(1,m);
x5=zeros(1,m);
for a=1:2
    if (a==1)
       x2(a)=x1(a)+x1(a+1);
       x3(a)=x1(a+2)+x1(a+3);
       x4(a)=x1(a+4)+x1(a+5);
       x5(a)=x1(a+6)+x1(a+7);
    else
       x2(a)=x1(a-1)-x1(a);
       x3(a)=x1(a+1)-x1(a+2);
       x4(a)=x1(a+3)-x1(a+4);
       x5(a)=x1(a+5)-x1(a+6);
    end
end
x6=zeros(1,n);
x7=zeros(1,n);
for b=1:2
    if (b==1)
        x6(b)=x2(b)+x3(b);
        x7(b)=x4(b)+x5(b);
    else
        x6(b)=x2(b)+((-1i).*x3(b));
        x7(b)=x4(b)+((-1i).*x5(b));
        
    end
end
for c=3:4
    if (c==3)
        x6(c)=x2(c-2)-x3(c-2);
        x7(c)=x4(c-2)-x5(c-2);
    else
        x6(c)=(x2(c-2)-((-1i).*x3(c-2)));
        x7(c)=(x4(c-2)-((-1i).*x5(c-2)));
    end
end
y=zeros(1,l);
for d=1:8
    g=exp(-1i*2*pi*(d-1)/8);
    h=exp(-1i*2*pi*(d-5)/8);
    if(d<=4)
        y(d)=x6(d)+(g.*(x7(d)));
    else
        y(d)=x6(d-4)-(h.*(x7(d-4)));
    end
end

clc;clear all;close all;
fd = 2000;
fs =1500;
Ts = 1/fs;
L=1000;
t = Ts*[0:L-1];
x = sin (2*pi*fd*t) 
fr = fs/L *(-L/2 : L/2 -1);
y=fft(x);
yabs=abs(fftshift(y));
fr1=fs/L*(-5*L/2:5*L/2-1);
yabs1=[yabs yabs yabs yabs yabs]
plot(fr1,yabs1)
xline(-fd,'--r')
xline(fd,'--r')

clc;clear all;close all;
a=[0:1:3];
x=[1 2 5 6];
b=[0:1:3];
c=-b;
h=[3 4 5 6]
m=length(x);
n=length(h);
X=[x,zeros(1,n)];
H=[h,zeros(1,m)];
n1=[-3:1:3];
Y=zeros(1,m+n-1);
for i=1:m+n-1
    for j=1:i
        if(j-i+m>0)
            Y(i)=Y(i)+X(j)*H(j-i+m);
        end
    end
end
figure();
subplot(2,2,1);stem(a,x,'filled');xlabel('n(samples)');ylabel('Amplitude(V)');
grid on;title('x[n]');
subplot(2,2,2);stem(b,h,'filled');xlabel('n(samples)');ylabel('Amplitude(V)');
grid on;title('h[n]');
subplot(2,2,3);stem(c,h,'filled');xlabel('n(samples)');ylabel('Amplitude(V)');
grid on;title('h[-n]');
subplot(2,2,4);stem(n1,Y,'filled');xlabel('n(samples)');ylabel('Amplitude(V)');
grid on;title('y[n]=x[n]*h[-n]');
sgtitle('Cross correlation  B.Mahaniya 2022105026')

clc;clear all;close all;
a=[0:1:3];
x=[1 3 5 7];
b=-a;
m=length(x);
X=[x,zeros(1,m)];
n=[-3:1:3];
Y=zeros(1,2*m-1);
for i=1:2*m-1
    for j=1:i
        if(j-i+m>0)
            Y(i)=Y(i)+X(j)*X(j-i+m);
        end
    end
end
figure();
subplot(3,1,1);stem(a,x,'filled');xlabel('n(samples)');ylabel('Amplitude(V)');
grid on;title('x[n]');
subplot(3,1,2);stem(b,x,'filled');xlabel('n(samples)');ylabel('Amplitude(V)');
grid on;title('x[-n]');
subplot(3,1,3);stem(n,Y,'filled');xlabel('n(samples)');ylabel('Amplitude(V)');
grid on;title('y[n]=x[n]*x[-n]');
sgtitle('Auto correlation  B.Mahaniya 2022105026')

clc;clear all; close all;
x=[1 3 5 7 2 4 6 8 9 1 2 3 4 5 6 7 1];
xle=length(x);
h=[1 2 3]
m= length(h);
L=5;
subplot(2,2,1)
stem(x);title('i/p x[n]')
subplot(2,2,2)
stem(h)

x=[x zeros(1,m-1)]
a=mod(length(x),L);
x=[x zeros(1,L-a)];
b=length(x)/L;
xb=reshape(x,L,b)'
xa=zeros(b,m-1)
for i=2:b
    for j=1:m-1
        xa(i,j)=xb(i-1,L-m+1+j)
    end
end

xb=[xa xb]
subplot(2,2,3)
stem(xb)

h=[h zeros(1,L-1)]
hk=fft(h);
subplot(3,1,2);
stem(hk);title('H(k)')
for i=1:b
    xk=fft(xb(i,:))
    yk(i,:)=xk.*hk
end
subplot(3,1,1);
stem(xk);title('X(k)')
subplot(3,1,3)
stem(yk);title('Y(k)')
for i=1:b
    yn(i,:)=ifft(yk(i,:))
end
y=[];
for i=1:b
    for j=m:(L+m-1)
    y=[y yn(i,j)]
    end 
end
y_n=[];
for i=1:xle+m-1
    y_n=[y_n y(i)]
end
figure
stem(uint8(y_n));
title('y(n)')

clc;clear all; close all;
x=[1 3 5 7 2 4 6 8 9 1 2 3 4 5 6 7 1];
xle=length(x);
h=[1 2 3]
m= length(h);
L=5;
subplot(2,2,1)
stem(x);title('i/p x[n]')
subplot(2,2,2)
stem(h)
x=[x zeros(1,m-1)]
a=mod(length(x),L)
x=[x zeros(1,L-a)]
b=length(x)/L;
xb=reshape(x,L,b)'
xa=zeros(b,m-1)
xb=[xb xa]
subplot(2,2,3)
stem(xb)

h=[h zeros(1,L-1)]
hk=fft(h);
figure
subplot(3,1,2)
stem(hk);title('H(k)')
for i=1:b
    xk=fft(xb(i,:))
    yk(i,:)=xk.*hk
end
subplot(3,1,1);
stem(xk);title('X(k)')
subplot(3,1,3)
stem(yk);title('Y(k)')
for i=1:b
    yn(i,:)=ifft(yk(i,:))
end
y=[];
for i=1:b
    for j=1:L+m-1
        if (j<=L)
            y = [y yn(i,j)]
        elseif(i<b)
            yn(i+1,j-L)=yn(i,j)+yn(i+1,j-L)
        end 
    end
end
y_n=[];
for i=1:xle+m-1
    y_n=[y_n y(i)]
end
figure
uint8(y_n)

clc;clear all;close all;
n=[-4:1:25]
a=n>=0;
b=n>=20;

A=7.*(a-b);subplot(2,2,1);xlabel('n(samples)');ylabel('Amplitude(V)');
stem(n,a,'filled');grid on;title('u[n]');
subplot(2,2,2);xlabel('n(samples)');ylabel('Amplitude(V)');
stem(n,b,'filled');grid on;title('u[n-20]');
subplot(2,2,3);xlabel('n(samples)');ylabel('Amplitude(V)');
stem(n,A,'filled');grid on;title('y[n]=7*(u[n]-u[n-20])');
sgtitle('Generation of Sequences');

clc;clear all;close all;
n=[0:10];
a=n>=0;
b=n>=4;
c=n.*a;
d=2.*b;
A=c-d;subplot(3,2,1);xlabel('n(samples)');ylabel('Amplitude(V)');
stem(n,a,'filled');grid on;title('u[n]');
subplot(3,2,2);xlabel('n(samples)');ylabel('Amplitude(V)');
stem(n,c,'filled');grid on;title('n*u[n]');
subplot(3,2,3);xlabel('n(samples)');ylabel('Amplitude(V)');
stem(n,b,'filled');grid on;title('u[n-4]');subplot(3,2,4);xlabel('n(samples)');ylabel('Amplitude(V)');
stem(n,d,'filled');grid on;title('2*u[n-4]');
subplot(3,2,5);xlabel('n(samples)');ylabel('Amplitude(V)');
stem(n,A,'filled');grid on;title('y[n]=n*u[n]-2*u[n-4]');
sgtitle('Generation of Sequences');

clc;clear all;close all;
n=[-10:1:10];
a=n==0;
b=n==4;
c=2.*b;
A=a-c;subplot(2,2,1);xlabel('n(samples)');ylabel('Amplitude(V)');
stem(n,a,'filled');grid on;title('i[n]');
subplot(2,2,2);xlabel('n(samples)');ylabel('Amplitude(V)');
stem(n,c,'filled');grid on;title('i[n-4]');
subplot(2,2,3);xlabel('n(samples)');ylabel('Amplitude(V)');
stem(n,c,'filled');grid on;title('2*i[n-4]');
subplot(2,2,4);xlabel('n(samples)');ylabel('Amplitude(V)');
stem(n,A,'filled');grid on;title('x[n]=i[n]-2*i[n-4]');
sgtitle('Generation of Sequences');

clc;clear all;close all;
n=[-10:1:10];
a=n==0;
b=n==4;
c=2.*b;
A=a-c;subplot(2,2,1);xlabel('n(samples)');ylabel('Amplitude(V)');
stem(n,a,'filled');grid on;title('i[n]');
subplot(2,2,2);xlabel('n(samples)');ylabel('Amplitude(V)');
stem(n,c,'filled');grid on;title('i[n-4]');
subplot(2,2,3);xlabel('n(samples)');ylabel('Amplitude(V)');
stem(n,c,'filled');grid on;title('2*i[n-4]');
subplot(2,2,4);xlabel('n(samples)');ylabel('Amplitude(V)');
stem(n,A,'filled');grid on;title('x[n]=i[n]-2*i[n-4]');
sgtitle('Generation of Sequences');

clc;clear all;close all;
x=[-2:0.01:5];
a=(4/3).*x;
b=x>=0;
c=x>=3;
d=b-c;
A=a.*d;subplot(3,2,1);plot(x,A);xlabel('time(s)');ylabel('Amplitude(V)');grid on;title('x(t)');
y1=x+3;subplot(3,2,2);plot(y1,A);xlabel('time(s)');ylabel('Amplitude(V)');grid on;title('x(t-3)');
y2=4.*x;subplot(3,2,3);plot(y2,A);xlabel('time(s)');ylabel('Amplitude(V)');grid on;title('x(t/4)');
y3=x/3;subplot(3,2,4);plot(y3,A);xlabel('time(s)');ylabel('Amplitude(V)');grid on;title('x(3t)');
y4=-x;subplot(3,2,5);plot(y4,A);xlabel('time(s)');ylabel('Amplitude(V)');grid on;title('x(-t)');

sgtitle('Continuous Signal');

stem(n,d,'filled');grid on;title('cos(0.12*pi*n)');

stem(n,d,'filled');grid on;title('cos(0.12*pi*n)');

A=a-b;subplot(2,2,1);stem(n,a,'filled');xlabel('n(samples)');ylabel('Amplitude(V)');
title('Generation of Sequences');

Interpolation and Decimation
clc;
clear all;
close all;
x = input('Enter the Input Sequence x(n): ');
I = input('Enter the Interpolation factor: ');
D = input('Enter the Decimation factor: ');
N = length(x);
n = 0:1:N-1;
subplot(311);
stem(n,x);
title('Input Sequence');
xlabel('n');
ylabel('Amplitude');

%Interpolation
y = [zeros(1,I*N)];
k=0;
for i = 1:1:I*N;
if mod(i-1,I)==0
k = k+1
y(i) = x(k);
else

2018105013
Benita .D
43

y(i) = 0;
end
end
j = 0:1:(I*N-1);
subplot(312);
stem(j,y);
title('Interpolated sequence');
xlabel('k');
ylabel('Amplitude');

%Decimation
p=0;
for i = 1:1:N
if mod(i-1,D)==0
p = p+1;
z(p) = x(i);
end
end
if mod((N/D),1)==0;
m = N/D-1;
else
m=floor(N/D);
end
n1=0:1:m;
subplot(313);
stem(n1,z);
title('Decimated sequence');
xlabel('k');
ylabel('Amplitude')

