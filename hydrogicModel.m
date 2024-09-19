%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Courtney DiVittorio and can be used to run a
% calibrated hydrologic model of the Sudd Wetland
% Last Updated Sept 19 2024 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all

%manually load variables that you would like to use into matlab
%The model was calibrated using TRMM precipitation and Hargreaves PET
%If using different variable names than what is described below modify
%accordingly

%% convert units to meters for all variables
P = P./1000; % precipitation [m] 
ET = PET./1000; % potential evapotranspiration [m]
A1obs = A1obs.*1000^2'; %[m^2]
Qin = Qin.*10^9; %[m^3]
Qout = Qout.*10^9; %[m^3/month]
Qoutmin = min(Qout(2:end)); 
A1min = min(A1obs(4:219));


%% set objective function weights
J = 10^10;
aJ1 = 1; %RMSE of area
aJ2 = 1; %RMSE of flow
aJ3 = 0.25; %timing of max flooded area
aJ4 = 0.25; %timing of min flooded area
aJ5 = 1; %bias in routing 


%% Simulate optimal - values provided in paper - supplemental data
clear S1sim S1simMid A1sim A1simMid Q12sim QoutSim S2sim S2simMid  
d = 1.123; %Optd;
OptrAind = 73;
rA = RechargeMatd02(:,OptrAind);
rR = 0;
aET = 0.945; 
aQout = 0.310; 
aQ12 = 0.1584; 
S1min = 2.5*10^9; 

%Initialize simulation variables
S1obs = S1min.*ones(180,1)+(A1obs(2:181)-A1min.*ones(180,1)).*d;
S1sim(1:180,1) = NaN;
S1sim(4,1) = S1obs(4); %m^3
S1sim(5,1) = S1obs(5);
A1sim(1:180,1) = NaN;
A1sim(4,1) = A1obs(4); %m^3
A1sim(5,1) = A1obs(4);
Q12sim(1:180,1) = NaN;
Q12sim(4,1)=aQ12*0.5*(S1sim(5)+S1sim(4));
QoutSim(1:180,1) = NaN;
QoutSim(4,1) = Qout(4);
%assume S2min
S2min = 0.25*S1min;
%1st value of S2 estimated from observed Qout
%and S-O relationship
S2sim(1:180,1) = NaN;
S2sim(4,1) = (1/aQout)*(0.5*(Qout(3)+Qout(4))-Qoutmin)+S2min;
S2sim(5,1) = (1/aQout)*(0.5*(Qout(4)+Qout(5))-Qoutmin)+S2min;
%to look at individual fluxes later
PmEAsim(1:180,1) = NaN;
RAsim(1:180,1) = NaN;
RRsim(1:180,1) = NaN;
%calculate beta and gamma terms
Beta0A = (1/(1+aQ12/2)).*(1+(1/(2*d)).*(P(1:180)-aET.*ET(1:180))-rA./d-aQ12/2); %Advance
Beta0R = (1/(1+aQ12/2)).*(1+(1/(2*d)).*(P(1:180)-aET.*ET(1:180))-rR/d-aQ12/2); %Recede
Beta1A = (1/(1+aQ12/2)).*((1/(2*d)).*(P(1:180)-aET.*ET(1:180))+rA./d); %Advance
Beta1R = (1/(1+aQ12/2)).*((1/(2*d)).*(P(1:180)-aET.*ET(1:180))+rR/d); %Recede
Beta2 = 1/(1+aQ12/2);
Beta3 = (1/(1+aQ12/2)).*(1/(2*d)).*(P(1:180)-aET.*ET(1:180)).*(2*d*A1min-2*S1min); 
%beta to plot
Beta3b = (1/(1+aQ12/2)).*(1/(2*d)).*(P(1:180)-aET.*ET(1:180)); 
Gamma0 = (1-aQout/2)/(1+aQout/2);
Gamma1 = 1/(1+aQout/2);
Gamma2 = (aQout*S2min-Qoutmin)/(1+aQout/2);
%simulate fluxes
Beta0 = NaN(180,1);
Beta1 = NaN(180,1);
for j = 6:180
   %simulate S1 - beg/end of month
   %need to compare previous 2 time steps b/c area 
   %is function of storage at prev time step
   if S1sim(j-1) > S1sim(j-2) %area advancing
       S1sim(j,1) = Beta0A(j-1)*S1sim(j-1)+Beta1A(j-1)*S1sim(j-2)...
           +Beta2*Qin(j-1)+Beta3(j-1);
       %store beta0 and beta1 so can plot
       Beta0(j-1,1) = Beta0A(j-1);
       Beta1(j-1,1) = Beta1A(j-1);
   else
       S1sim(j,1) = Beta0R(j-1)*S1sim(j-1)+Beta1R(j-1)*S1sim(j-2)...
           +Beta2*Qin(j-1)+Beta3(j-1);
       Beta0(j-1,1) = Beta0R(j-1);
       Beta1(j-1,1) = Beta1R(j-1);
   end
   %simulate Q12 and A1 for time j
   Q12sim(j-1) = aQ12*0.5*(S1sim(j)+S1sim(j-1)); %middle of month
   A1sim(j) = (1/d)*(S1sim(j-1)-S1min+d*A1min); %beg/end of month
   %simulate S2 
   S2sim(j,1)=Gamma0*S2sim(j-1)+Gamma1*Q12sim(j-1)+Gamma2; %Beg/end of month
   %solve for Qout
   QoutSim(j-1)=Qoutmin+aQout*(0.5*(S2sim(j)+S2sim(j-1))-S2min); %middle of month
   %A(P-E) to plot later
   PmEAsim(j-1,1) =(P(j-1)-aET*ET(j-1))*(0.5*(A1sim(j)+A1sim(j-1))); %middle of month
   %recharge to plot later - need to get rid of
   %negative values at end
   RAsim(j-1,1)=rA(j-1)*(A1sim(j)-A1sim(j-1)); %middle of month
   RRsim(j-1,1)=rR*(A1sim(j-1)-A1sim(j)); %middle of month
end
%recalculate S2min and iterate until convergence
%affects S2 simulation and Qout, measure in
%terms of change in Qout since this is what we
%are calibrating to
RAsim(RAsim<0)=NaN;
RRsim(RRsim<0)=NaN;
percDiffQ = 1;
kp=0;
while (percDiffQ > 0.001) && kp < 31
    kp = kp+1;
    S2minPrev = S2min;
    S2min = min(S2sim);
    %new initial value S2
    S2sim(4,1) = (1/aQout)*(0.5*(Qout(3)+Qout(4))-Qoutmin)+S2min;
    S2sim(5,1) = (1/aQout)*(0.5*(Qout(4)+Qout(5))-Qoutmin)+S2min;
    %recalc gamma2
    Gamma2 = (aQout*S2min-Qoutmin)/(1+aQout/2);
    QoutPrev = QoutSim;
    for j = 6:180
        %simulate S2 routing component
        S2sim(j,1)=Gamma0*S2sim(j-1)+Gamma1*Q12sim(j-1)+Gamma2;
        %solve for Qout
        QoutSim(j-1)=Qoutmin+aQout*(0.5*(S2sim(j)+S2sim(j-1))-S2min);
    end
    %check change in Q at each time step, cant
    %exceed 0.1%
    percDiffQ = max(abs(QoutSim(5:179)-QoutPrev(5:179))./QoutSim(5:179));
    %diplay if kp reaches 30
    if kp == 30 && percDiffQ > 0.001
        disp('S2min did not converge')
    end
end
%interpolate and find flooded area at
%middle of month to compare to observed
S1simMid(:,1)= 0.5.*(S1sim(1:end-1)+S1sim(2:end));
A1simMid(:,1) = 0.5.*(A1sim(1:end-1)+A1sim(2:end));
S2simMid(:,1)=0.5.*(S2sim(1:end-1)+S2sim(2:end));

A1simMid(:,1) = 0.5.*(A1sim(1:end-1)+A1sim(2:end));
 %evaluate objective function
%RMSE
Jtmp(1,1) =  aJ1*((mean(((A1obs(5:173)-A1simMid(5:173))/...
   (max(A1obs)-min(A1obs))).^2))^0.5);
Jtmp(2,1) =  aJ2*((mean(((Qout(5:173)-QoutSim(5:173))/...
   (max(Qout)-min(Qout))).^2))^0.5);
%difference in timing between max and min flooded area and
max1a = int8(ones(14,1));
min1a = int8(ones(14,1));
max1b = int8(ones(14,1));
min1b  = int8(ones(14,1));
for yr = 1:14
    tmp1a = A1obs(12*(yr-1)+6:12*yr+6); %observed
    tmp1b =  A1simMid(12*(yr-1)+6:12*yr+6); %simulated
    tmpmax = find(tmp1a == max(tmp1a));
    max1a(yr) = tmpmax(1); 
    tmpmin = find(tmp1a == min(tmp1a));
    min1a(yr) = tmpmin(1); 
    tmpmax = find(tmp1b == max(tmp1b));
    max1b(yr) = tmpmax(1); 
    tmpmin = find(tmp1b == min(tmp1b));
    min1b(yr) = tmpmin(1);     
end
%for min values replace 1 with 13 and 2 with 14 so
%can add and subtract accross multiple years
min1a(min1a == 1) = 13;
min1b(min1b == 1) = 13;
min1a(min1a == 2) = 14;
min1b(min1b == 2) = 14;
Jtmp(3,1) = aJ3*((mean(((max1a-max1b)/2).^2))^0.5);
Jtmp(4,1) = aJ4*((mean(((min1a-min1b)/4).^2))^0.5);
%bias between routing and outflows over full
%simulation
Jtmp(5,1) = aJ5*(abs(sum(Q12sim(6:179))-sum(Qout(6:179)))...
   /(0.05*sum(Qin(6:179)-Qout(6:179))));
%record optimal A and Q sequences
%OptA1sim_calQA = A1simMid;
%OptQoutSim_calQA = QoutSim;

% OptA1sim_calQ = A1simMid;
% OptQoutSim_calQ = QoutSim;

%% record flow and area simulation for best estimate
QoutSim1 = QoutSim;
A1simMid1 = A1simMid;
S1simMid1 = S1simMid;
S2simMid1 = S2simMid;
%NSE = 0.7665
%% record Asim for a only simulation
A1simMid2 = A1simMid;

%% record Qoutsim for q only
QoutSim2 = QoutSim;
%NSE = 0.8276
%% Plots
Dates = SuddFlowEstimates.Date;
j=0;
for yr = 2001:2014
    tmpInd = find(year(Dates) == yr & month(Dates) == 1);
    j=j+1;
    DateTick(j,1) = (Dates(tmpInd));
    DateLabels{j} = year(Dates(tmpInd));
end

%stats and bias
NSE = 1-sum((QoutSim(6:179) - Qout(6:179)).^2)/sum((Qout(6:179)-mean(Qout(6:179))).^2);
%
sum(Q12sim(6:173))
sum(Qout(6:173))
sum(QoutSim(6:173))
NSEa = 1-sum((A1simMid(5:173) - A1obs(5:173)).^2)/sum((A1obs(5:173)-mean(A1obs(5:173))).^2);

%% areas only
figure
plot(Dates(4:180-7),A1obs(4:180-7)./1000^2,...
    'LineWidth',1.5,'Color',[0 76 153]./255)
hold on
plot(Dates(4:180-7),A1simMid(4:180-7)./1000^2,...
    'LineWidth',1.5,'Color',[161 25 25]./255)

ax = gca;
ax.YRuler.Exponent = 0;
ax.FontSize = 12;
ytickformat('%,4.0g')
lg = legend('Estimated','Simulated');
lg.FontSize = 11;
title('Comparison between Flooded Areas for Smaller Depth','FontSize',14)
grid on
ylabel('Flooded Area (km^{2})','FontSize',13)
xticks(DateTick)
xticklabels(DateLabels)
%% plot simulated and observed areas and outflows
figure('Position',[100 -300 1200 800])
subplot(3,1,1)
plot(Dates(4:180-7),A1obs(4:180-7)./1000^2,...
    'LineWidth',1.5,'Color',[0 76 153]./255)
hold on
plot(Dates(4:180-7),A1simMid1(4:180-7)./1000^2,...
    'LineWidth',1.5,'Color',[161 25 25]./255)
hold on
plot(Dates(4:180-7),A1simMid2(4:180-7)./1000^2,...
    'LineWidth',1.5,'Color',[161 25 25]./255,'LineStyle','--')
ax = gca;
ax.YRuler.Exponent = 0;
ax.FontSize = 12;
ytickformat('%,4.0g')
lg = legend('Estimated - Satellite','Simulated - Joint Calibration',...
    'Simulated - Area Only Calibration','location','southoutside','orientation','horizontal');
lg.FontSize = 11;
title('Comparison between Flooded Areas','FontSize',14)
grid on
ylabel('Flooded Area (km^{2})','FontSize',13)
xticks(DateTick)
xticklabels(DateLabels)
%ylim([12000 27000])
% hold on
% yyaxis right
% plot(Dates(4:180-7),P(4:180-7).*1000,...
%     'LineWidth',1.5)

subplot(3,1,2)
plot(Dates(4:180-7),Qout(4:180-7)./10^9,...
    'LineWidth',1.5,'Color',[0 76 153]./255)
hold on
plot(Dates(4:180-7),QoutSim1(4:180-7)./10^9,...
    'LineWidth',1.5,'Color',[161 25 25]./255)
hold on
plot(Dates(4:180-7),QoutSim2(4:180-7)./10^9,...
    'LineWidth',1.5,'Color',[161 25 25]./255,'LineStyle','--')
set(gca,'FontSize',12)
lg = legend('Estimated - Regression Models','Simulated - Joint Calibration',...
    'Simulated - Flow Only Calibration','location','southoutside','orientation','horizontal');
lg.FontSize = 11;
title('Comparison between Outflows','FontSize',14)
grid on
ylabel('Outflow (bcm/month)','FontSize',13)
xticks(DateTick)
xticklabels(DateLabels)

subplot(3,1,3)
plot(Dates(4:180-7),S1simMid(4:180-7)./10^9,...
    'LineWidth',1.5,'Color',[6 191 167]./255)
set(gca,'FontSize',12)
grid on
ylabel('Storage (bcm/month)','FontSize',13)
xticks(DateTick)
xticklabels(DateLabels)
hold on
plot(Dates(4:180-7),S2simMid(4:180-7)./10^9,...
    'LineWidth',1.5,'Color',[86 30 156]./255)
grid on
ylabel('Storage (bcm/month)','FontSize',13)
xticks(DateTick)
xticklabels(DateLabels)
lg = legend('S1 Storage - Joint Calibration',...
    'S2 Storage - Joint Calibration','location','southoutside','orientation','horizontal');
lg.FontSize = 11;
set(gca,'FontSize',12)
title('Simulated Storages','FontSize',16)





%% scatter plots
figure %('Position',[100 20 1150 1150])
subplot(2,2,1)
scatter(A1simMid(4:173)./1000^2,A1obs(4:173)./1000^2)
minBound = min([min(A1simMid)/1000^2-0.05*(max(A1simMid)/1000^2-min(A1simMid)/1000^2),...
    min(A1obs)-0.05*(max(A1obs)/1000^2-min(A1obs)/1000^2)]);
maxBound = max([max(A1simMid)/1000^2+0.05*(max(A1simMid)/1000^2-min(A1simMid)/1000^2),...
    max(A1obs)/1000^2+0.05*(max(A1obs)/1000^2-min(A1obs)/1000^2)]); 
xlim([minBound, maxBound])
ylim([minBound, maxBound])
hold on
plot([minBound,maxBound],[minBound,maxBound],'LineWidth',1.25,...
    'Color',[205 48 48]./255,'LineStyle','--')
grid on
ytickformat('%,4.0g')
xtickformat('%,4.0g')
ax = gca;
ax.YRuler.Exponent = 0;
ax.XRuler.Exponent = 0;
ax.FontSize = 12;
xlabel('Simulated Flooded Area (km^{2})','FontSize',13)
ylabel('Estimated Flooded Area (km^{2})','FontSize',13)
title('Scatter Plot of Flooded Areas','FontSize',15)
lg = legend('Data','1-1 Line');
lg.FontSize = 12;
box on

subplot(2,2,2)
scatter(QoutSim(4:173)./10^9,Qout(4:173)./10^9)
minBound = min([min(QoutSim)/10^9-0.05*(max(QoutSim)/10^9-min(QoutSim)/10^9),...
    min(Qout)-0.05*(max(Qout)/10^9-min(Qout)/10^9)]);
maxBound = max([max(QoutSim)/10^9+0.05*(max(QoutSim)/10^9-min(QoutSim)/10^9),...
    max(Qout)/10^9+0.05*(max(Qout)/10^9-min(Qout)/10^9)]); 
xlim([minBound, maxBound])
ylim([minBound, maxBound])
hold on
plot([minBound,maxBound],[minBound,maxBound],'LineWidth',1.25,...
    'Color',[205 48 48]./255,'LineStyle','--')
grid on
set(gca,'FontSize',12)
xlabel('Simulated Outflow (bcm)','FontSize',13)
ylabel('Estimated Outflow (bcm)','FontSize',13)
title('Scatter Plot of Outflows','FontSize',15)
lg = legend('Data','1-1 Line');
lg.FontSize = 12;
box on
%% plot fitted relaitonships for storage-area and outflow-storage
%figure('Position',[100 -300 1400 700])
subplot(2,2,3)
scatter((A1obs(2:174)-A1min)./10^6,(S1simMid(1:173))./10^9)
xv = [0 nanmax(A1obs)-nanmin(A1obs)].*ones(1,2);
yv = S1min.*ones(1,2)+d.*(xv);
hold on
plot(xv./10^6,yv./10^9,'r')
xtickformat('%,4.0g')
ax = gca;
ax.XRuler.Exponent = 0;
ax.FontSize = 12;
title('Relationship Between S1 Storage and Area','FontSize',15)
grid on
xlabel('A_{OBS}(k) km^{2}','FontSize',13)
ylabel('S1(k-1) bcm','FontSize',13)
lg = legend('Data','Fitted Linear Relationship');
lg.FontSize = 12;
box on

subplot(2,2,4)
scatter((S2simMid(1:173)- S2min)./10^9,(Qout(1:173))./10^9)
xv = [0 nanmax(S2simMid)-nanmin(S2simMid)];
yv = Qoutmin.*ones(1,2)+aQout.*(xv);
hold on
plot(xv./10^9,yv./10^9,'r')
set(gca,'FontSize',12)
title('Relationship Between Outflow and S2 Storage','FontSize',15)
grid on
xlabel('S2(k) bcm','FontSize',13)
ylabel('Q_{OUT,OBS}(k) bcm','FontSize',13)
lg = legend('Data','Fitted Linear Relationship');
lg.FontSize = 12;
box on

%% For A only
figure('Position',[100 -300 1400 700])
subplot(1,2,1)
scatter(A1simMid(4:173)./1000^2,A1obs(4:173)./1000^2)
ytickformat('%,4.0g')
xtickformat('%,4.0g')
ax = gca;
ax.YRuler.Exponent = 0;
ax.XRuler.Exponent = 0;
ax.FontSize = 12;
minBound = min([min(A1simMid)/1000^2-0.05*(max(A1simMid)/1000^2-min(A1simMid)/1000^2),...
    min(A1obs)-0.05*(max(A1obs)/1000^2-min(A1obs)/1000^2)]);
maxBound = max([max(A1simMid)/1000^2+0.05*(max(A1simMid)/1000^2-min(A1simMid)/1000^2),...
    max(A1obs)/1000^2+0.05*(max(A1obs)/1000^2-min(A1obs)/1000^2)]); 
xlim([minBound, maxBound])
ylim([minBound, maxBound])
hold on
plot([minBound,maxBound],[minBound,maxBound],'LineWidth',1.25,...
    'Color',[205 48 48]./255,'LineStyle','--')
grid on
xlabel('Simulated Flooded Area (km^{2})','FontSize',13)
ylabel('MODIS-derived Flooded Area (km^{2})','FontSize',13)


title('Scatter Plot of Flooded Areas','FontSize',16)
lg = legend('Data Points','1-1 Line');
lg.FontSize = 12;
box on

% plot fitted relaitonships for storage-area and outflow-storage
subplot(1,2,2)
scatter((A1obs(2:174)-A1min)./10^6,(S1simMid(1:173))./10^9)
xv = [0 nanmax(A1obs)-nanmin(A1obs)].*ones(1,2);
yv = S1min.*ones(1,2)+d.*(xv);
hold on
plot(xv./10^6,yv./10^9,'r','LineWidth',1.25,...
    'Color',[205 48 48]./255,'LineStyle','--')
grid on
xlabel('A_{OBS}(k) - A_{MIN} (km^{2})','FontSize',13)
xtickformat('%,4.0g')
ax = gca;
ax.FontSize = 12;
ax.XRuler.Exponent = 0;
ylabel('S1(k-1) (bcm/month)','FontSize',13)
lg = legend('Data Points','Fitted Linear Relationship');
lg.FontSize = 12;
box on
title('Relationship Between S1 Storage and Area','FontSize',16)

%% beta's and gammas

figure('Position',[100 -300 1200 800])
subplot(2,1,1)
hold all
plot(Dates(4:180-7),Beta0(4:180-7),...
    'LineWidth',1.5,'Color',[15 135 148]./255)
plot(Dates(4:180-7),Beta1(4:180-7),...
    'LineWidth',1.5,'Color',[158 53 8]./255)
plot(Dates(4:180-7),Beta2.*ones(170,1),...
    'LineWidth',1.5,'Color',[27 117 31]./255)
plot(Dates(4:180-7),Beta3b(4:180-7),...
    'LineWidth',1.5,'Color',[114 9 133]./255)
ax = gca;
%ax.YRuler.Exponent = 0;
ax.FontSize = 12;
%ytickformat('%,4.0g')
lg = legend('\beta_{0}','\beta_{1}','\beta_{2}','\beta_{3}');
lg.FontSize = 12;
lg.Orientation = ('horizontal');
title('\beta Terms for S1 Simulation','FontSize',16)
grid on
ylabel('\beta Value','FontSize',13)
xticks(DateTick)
xticklabels(DateLabels)
%ylim([12000 27000])
box on


