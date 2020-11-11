%% Reward Sensitive Gaussian Process 
% "Reinforcement regulates timing variability" 
% Jing Wang  jingwang.physics(a)gmail.com
%  - March 2019



%%  Parameters
clc; clear; close all;
       Nmax = 10e3;
    Wp      = .1;  % weber constant
    sig_g   = 2*Wp; % total variance
    l1      = 20; % long range correlation length, in unit of trials
    l2      = 2; % reward sensitive correlation length, in unit of trials
    Alpha   = [1.0, 0.0, 0.5]; % fraction of variance for slow component,  0 ~ 1
    NAlpha  = length(Alpha);
    sig_g0  = 1*Wp;  % variance of private noise
    acLength= ceil(2.5*max(l1,l2)); % max length of correlation to be considered
    k_reward= 0.1; % reward acceptance window
   dk_reward= 0.01;
k_reward_min= 0.05; 
k_reward_max= 0.3;
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
disp('++++++++++++++++ RSGP simulation +++++++++++++++++'); 

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

% Mask matrix for reward
  MK = ones(length(Reward),length(Reward)); 
inde = find(Reward ==0);   % error trials
        for k = 1:length(inde)
                  kd = inde(k)+1; %min(indx(k)+1,size(K,1));
                  if kd<size(MK,2)
                   MK(inde(k),kd:end) = 0;  MK(kd:end,inde(k)) = 0; 
                  end

        end
 MK = (MK+MK')./2 ; 

%% Maxmium Marginal Likelihood fit to model parameters

disp('++++++++++++ RSGP Hyperparameters Fit ++++++++++++');
disp(' *Warning* Turtle ahead - Continue with fitting?');
prompt=' Y/N ?';

FIT = input(prompt,'s');

if FIT =='y' | FIT =='Y'
    randn('state',0);
    X = [1:Nmax]';
    Nsample = 50;
    
    % step one: slow
    meanfunc =  {@meanConst} ;                   %  don't use a mean function
    covfunc  = {@covSum, {'covSEiso'}};% feval(covfunc{:}), covariance function
    likfunc  = @likGauss;   % Gaussian likelihood
    inf = @infGaussLik;
    
    true_para1=[l1; sig_g.*sqrt(alpha)];
    hyp = struct('mean', [0], 'cov', log(true_para1), 'lik', -1);
    ind =(Reward(1:end-1)==0);
    hyp_opt1 = minimize(hyp, @gp, -1e3, inf, meanfunc, covfunc, likfunc, X(ind), Tp(ind)); %  maxlikeli
    hyp.mean =  hyp_opt1.mean;
    hyp.lik =  hyp_opt1.lik;
    
    test_para1=[];
    for i=1:length(true_para1) % grid
        test_para1(:,i) = true_para1(i) + 2*true_para1(i).*[-ceil(Nsample/2)+1:1:floor(Nsample/2)]'./Nsample;
    end
    D = size(test_para1,1);
    negLL1 = nan(D,D);
    for i=1:D
        for j=1:D
            disp(num2str([i,j] ));
            hyp.cov = log([test_para1(i,1);test_para1(j,2)]);
            negLL1(i,j) = gp(hyp, inf, meanfunc, covfunc, likfunc, X(ind), Tp(ind));
        end
    end
    
    % step two: slow and fast
    meanfunc =  {@meanConst} ;
    covfunc  = {@covSum, {'covSEiso','covSEisoRew','covNoise'}};%{@covSum, {'covSEiso','covSEisoRew','covNoise'}}; feval(covfunc{:}), covariance function
    likfunc  = @likGauss;
    
    prior.mean={[]}; % with prior
    prior.cov ={{'priorDelta'};{'priorDelta'};[];[];[]};
    prior.lik ={[]};
    inf = {@infPrior, @infGaussLik, prior};
    true_para2 =[l2; sig_g.*sqrt(1-alpha)];
    init_cov=[exp(hyp_opt1.cov); true_para2 ;sig_g0 ]; % initial parameters, use result of previous step
    hyp_init2 = struct('mean', hyp_opt1.mean, 'cov', log(init_cov), 'lik', hyp_opt1.lik);
    
    global MaskK
    MaskK = MK;
    hyp_opt2 = minimize(hyp_init2, @gp, -1e3, inf, meanfunc, covfunc, likfunc, X, Tp); %  maxlikeli
    hyp.mean =  hyp_opt2.mean;
    hyp.lik =  hyp_opt2.lik;
    
    test_para2=[];
    for i=1:(length(true_para2)) % grid
        test_para2(:,i) = true_para2(i) + true_para2(i).*[-ceil(Nsample/2)+1:1:floor(Nsample/2)]'./Nsample;
    end
    D = size(test_para2,1);
    negLL2 = nan(D,D);
    for i=1:D
        for j=1:D
            disp(num2str([i,j] ));
            hyp.cov = [hyp_opt2.cov(1:2); log([test_para2(i,1);test_para2(j,2)]);hyp_opt2.cov(end)];
            negLL2(i,j) = gp(hyp, inf, meanfunc, covfunc, likfunc, X, Tp);
        end
    end
    
    figure('Name','Likelihood profile','Position',[ 0 0 450 900]);
    subplot(211)
    imagesc(test_para1(:,1),test_para1(:,2),(-negLL1)); hold on
    plot(true_para1(1),true_para1(2),'b*','MarkerSize',20); hold on
    %plot(exp(hyp_opt1.cov(1)),exp(hyp_opt1.cov(2)),'k*','MarkerSize',10); hold on
    xlabel('l1');ylabel('\sigma1');drawnow;
    subplot(212)
    imagesc(test_para2(:,1),test_para2(:,2),(-negLL2)); hold on
    plot(true_para2(1),true_para2(2),'b*','MarkerSize',20); hold on
    %plot(exp(hyp_opt2.cov(3)),exp(hyp_opt2.cov(4)),'k*','MarkerSize',10); hold on
    xlabel('l2');ylabel('\sigma2');drawnow;
end
