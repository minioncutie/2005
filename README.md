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


