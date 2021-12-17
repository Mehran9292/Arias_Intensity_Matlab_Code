clc
clear
ACC=importdata('AGW000_drift.txt');
ACC=ACC.data;
Q = permute(ACC,[2 1]);
A=Q(:);
Acceleration=A*9.81;
deltaT=0:0.005:39.995;
m=length(Acceleration);
z=zeros(m,2);
z(:,1)=deltaT(1,:);
z(:,2)=Acceleration(:,1);
z(1,1)=0;
z(1,2)=0;
UncorrectedA=z(:,2);
t=z(:,1);
Uncorrectedvelocity = cumtrapz(deltaT,z(:,2));
figure(1)
subplot(3,1,1)
%plot(z(:,1),UncorrectedA/9.81)
hold on
plot(z(:,1),CorrectedA/9.81)
xlabel('Time (sec)') 
ylabel('Acceleration (m/s^2)') 
legend('Uncorrected','Corrected')
subplot(3,1,2)
%plot(z(:,1),Uncorrectedvelocity(:,1)*100)
hold on
plot(z(:,1),Correctedvelocity(:,1)*100)
xlabel('Time (sec)') 
ylabel('Velocity (cm/s)')
legend('Uncorrected','Corrected')
subplot(3,1,3)
%plot(z(:,1),Uncorrecteddisplacement(:,1)*100)
hold on
plot(z(:,1),CorrectedDisplacement(:,1)*100)
xlabel('Time (sec)') 
ylabel('Displacement (cm)')
legend('Uncorrected','Corrected')
% T=(0.01:0.01:3)';
% figure(2)
% G = 9.80665;
%Intensity = cumsum(UncorrectedA.^2)*pi*0.05/2/G;
% Intensity = cumtrapz(UncorrectedA.^2)*pi*0.05;
% Intensity2 = cumsum(UncorrectedA.^2)*pi*40;
% Ht=Intensity./Intensity2;
% plot(z(:,1), Ht)

f=0;
for i=1:length(deltaT)-1 
    f=f+((((UncorrectedA(i+1,1).^2)+(UncorrectedA(i,1).^2)).*0.005)./2);
end 

%xs=sum(H);
f2=0;
for i=1:length(deltaT)-1 
    f2(i+1)=f2(i)+(((((UncorrectedA(i+1,1).^2)+(UncorrectedA(i,1).^2)).*0.005)./2)./f);
end 
figure(2)
plot(z(:,1),f2)


function ai = arias_intensity(t,UncorrectedA)
G = 9.80665;
ai = trapz(t, (UncorrectedA.*UncorrectedA));
ai = ai * (pi / (2 * G));
end

