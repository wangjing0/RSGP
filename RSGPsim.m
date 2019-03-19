%% RSGP simulating

%%  Parameters
clc; clear; close all;
       Nmax = 1e4;
    Wp      = .1;
    sig_g   = 2*Wp;
    l1      = 20; % correlation length, in unit of trials
    l2      = 1;
    Alpha   = [1, 0.5, 0.0];
    NAlpha  = length(Alpha);
    sig_g0  = 1*Wp;  
    acLength= ceil(2.5*max(l1,l2)); % max length of correlation to be considered
    k_reward= 0.1; 
   dk_reward= 0.01;
k_reward_min= 0.02; 
k_reward_max= 0.5;
 kernelfunc = 'squaredexponential';

     N_bins = 8;
       blow = 1; 
       bhigh= 5;
        Bin = [-1*2.^(linspace(-blow,-bhigh,ceil(N_bins/2))),0,fliplr(2.^(linspace(-blow,-bhigh,ceil(N_bins/2))))];
   maxNlags = 3*max(l1,l2);  Nlags = 1:maxNlags;

%% GP kernel
[covInd1,covInd2] = meshgrid((1:Nmax),(1:Nmax));

switch kernelfunc
    case 'squaredexponential'
        cov_g1 = sig_g^2 * exp(- (covInd1-covInd2).^2 ./(2.*l1.^2)) ;
        cov_g2 = sig_g^2 * exp(- (covInd1-covInd2).^2 ./(2.*l2.^2)) ;
        cov_g3 = sig_g0^2.* eye(Nmax);
    case 'exp'
        cov_g = sig_g^2 * exp(- (abs(covInd1-covInd2)) ./l1) + sig_g0^2.* eye(Nmax);
    case  'equal'
        cov_g = sig_g^2 * ones(size(covInd1)) + sig_g0^2.* eye(Nmax);
    case 'whitenoise'
        cov_g = sig_g^2.* eye(Nmax);
    case 'besselk'
        nu=l-1/2; % order of, modified Bessel function of second kind
        xi_xj = sqrt(2.*nu).*abs(covInd1-covInd2)./l1;
        cov_g = sig_g^2.*(1./gamma(nu)./2.^(nu-1)).*((xi_xj).^nu).*besselk(nu,xi_xj);
        cov_g(logical(eye(Nmax)))=sig_g^2 + sig_g0^2;
end


    
%% Simulation and plotting
  

h1 = figure('Name','Tp','Position',[ 0 500 1200 900]);
h2 = figure('Name','MuSigma(Tp[n]) - Tp[n-1]','Position',[ 0 500 1200 300]);
h3 = figure('Name','ParCorr','Position',[ 0 0 1200 500]);
h4 = figure('Name','cov','Position',[ 0 0 1200 300]);
for iii=1:NAlpha
    alpha = Alpha(iii);    
 
    Tp=nan(Nmax,1);
    mTp=nan(Nmax,1); mTp(1)= 0;
    sigTp=nan(Nmax,1); sigTp(1)=sig_g^2;
    Reward=nan(Nmax,1);
    randn('state',1);
    
    rand1=randn(Nmax,1);
    rand2=randn(Nmax,1);
    
    for i=1:Nmax
        Tp(i) = mTp(i) +  rand1(i).*sqrt(sigTp(i) ) ;
    
        Reward(i) = abs(Tp(i)) < k_reward ;%(abs(Tp(i)-TS)< TS*k_reward ).* (1 - abs(Tp(i)-TS)./(TS*k_reward)); % analog reward
        
        if dk_reward % reward window on staircase
            if Reward(i)>0
                k_reward=max(k_reward_min, k_reward-dk_reward);
            else
                k_reward=min(k_reward_max, k_reward+dk_reward);
            end
        end
        
        if i < Nmax
            clear K* Tp_
            K1  = cov_g1(max(1,i-acLength):i,max(1,i-acLength):i); % ixi cov of known
            K_1 = cov_g1(i+1,max(1,i-acLength):i); % 1xi
            K__1= cov_g1(i+1,i+1);  % 1x1 cov of predict
            K2  = cov_g2(max(1,i-acLength):i,max(1,i-acLength):i); % ixi cov of known
            K_2 = cov_g2(i+1,max(1,i-acLength):i); % 1xi
            K__2= cov_g2(i+1,i+1);  % 1x1 cov of predict
            K3  = cov_g3(max(1,i-acLength):i,max(1,i-acLength):i);
            K__3= cov_g3(i+1,i+1);
            Tp_ = Tp(max(1,i-acLength):i);
            
                indx= find(Reward(max(1,i-acLength):i)==0);% only use rewarded trials for update
              for k = 1:length(indx)
                  kd = min(indx(k)+1,size(K2,1));
                K2(indx(k),kd:end) = 0;  K2(kd:end,indx(k)) = 0; 
              end
                K2(logical(eye(size(K2))))=cov_g2(i,i) ;
                K_2(indx) = 0;
                
                 sigTp(i+1)= (alpha.*K__1+ (1-alpha).*K__2 + K__3) - ...
                     (alpha*K_1+ (1-alpha).*K_2)*inv(alpha*K1 + (1-alpha).*K2 + K3)*(alpha*K_1 + (1-alpha).*K_2)' ;

                 mTp(i+1)  =real( (alpha.*K_1 + (1-alpha).*K_2)*inv(alpha*K1 + (1-alpha).*K2 + K3) * Tp_) ;                
       
        end        
    end

        
    time = [1:Nmax]'; npts = min(1e3,Nmax);
    indic = boolean(zeros(size(Tp))); indic(Nmax-npts+1:Nmax)=true;
    figure(h1);
    sh1(iii)=subplot(NAlpha,1,iii);
    shadedErrorBar(time(indic),mTp(indic),sqrt(sigTp(indic)),'lineprops',{'k-','LineWidth',1,'markerfacecolor','k'}); hold on;drawnow
    plot(time(Reward==0 &indic ),Tp(Reward==0 & indic ),'k.'); hold on
    plot(time(Reward>0 & indic ),Tp(Reward>0 & indic ),'g.');drawnow;makeaxis();

    z = Tp;
    z_1 = [nan;z(1:end-1)];
    
    mTp_stdTp=nan(length(Bin)-1,3);
    [pc,~,cf] = parcorr(z,maxNlags,[],2);
    for k = 1:length(Bin)-1
        ind = (z_1>= Bin(k) &  z_1 < Bin(k+1));
        if sum(ind)
            mTp_stdTp(k,:)=[nanmean(z_1(ind)),nanmean(z(ind)),nanstd(z(ind))];
        end
    end
    
    figure(h2);sh2(iii)=subplot(1,NAlpha,iii);
    plot(mTp_stdTp(:,1),mTp_stdTp(:,3),'b-','LineWidth',2); hold on;
    plot(mTp_stdTp(:,1),mTp_stdTp(:,3),'wo','MarkerSize',12,'MarkerFaceColor','b'); hold on;
    plot(mTp_stdTp(:,1),mTp_stdTp(:,2),'k-','LineWidth',2); hold on;
    plot(mTp_stdTp(:,1),mTp_stdTp(:,2),'wo','MarkerSize',12,'MarkerFaceColor','k'); hold on;
    plot(mTp_stdTp(:,1),zeros(size(mTp_stdTp(:,1))),'k:');
    drawnow;makeaxis();
    
    figure(h3);sh3(iii)=subplot(1,NAlpha,iii);
    plot(Nlags,pc(2:end),'k-','LineWidth',2); hold on;
    plot(Nlags,ones(size(Nlags)).*cf,'k:');hold on;drawnow; %makeaxis();

end

figure(h4);
for jj=1:NAlpha
    sh4(jj)=subplot(1,NAlpha,jj);
    imagesc(log((Alpha(jj)*K1 + (1-Alpha(jj)).*K2 + K3)));axis equal;axis off
end
figure(h1);linkaxes(sh1,'xy');
figure(h2);linkaxes(sh2,'xy');
figure(h3);linkaxes(sh3,'xy');
figure(h4);linkaxes(sh4,'xy');