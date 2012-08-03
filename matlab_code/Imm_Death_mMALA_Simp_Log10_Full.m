function [ ] = Imm_Death_mMALA_Simp_Log10_Full( y, Time, k1, k2 )

warning off

% Random Numbers...
randn('state', sum(100*clock));
rand('twister', sum(100*clock));

[N,D] = size(y);
D = 2;
        
TimeStep = diff(Time);

% Langevin Setup
NumOfIterations = 10000;
BurnIn          = 5000;
StepSize        = 0.05;


Proposed = 0;
Accepted = 0;


% Set initial values of w
Paras  = [log10(k1) log10(k2)];
ParasSaved = zeros(NumOfIterations-BurnIn,D);

% Calculate joint log likelihood for current w
%LogPrior      = LogNormPDF(zeros(1,D),Paras,alpha);

Mu  = Calc_Mean(k1,k2,y(1:N-1),TimeStep);
Var = Calc_Variance(k1,k2,y(1:N-1),TimeStep);
%LogLikelihood  = sum( log( normpdf(y(2:end),Mu,sqrt(Var)) ) );
LogLikelihood = pn(10^k1, 10^k2, TimeStep, y);

CurrentLJL    = LogLikelihood; % + LogPrior;


%%% Pre-Calculations %%%
ParasNew = Paras;

% Calculate G based on new parameters
% Calculate metric tensor
dM_dk1     = dMean_dk1(k1,k2,y(1:N-1),TimeStep);
dSigma_dk1 = dVariance_dk1(k1,k2,y(1:N-1),TimeStep);
dSigma_dk1 = diag(dSigma_dk1);

dM_dk2     = dMean_dk2(k1,k2,y(1:N-1),TimeStep);
dSigma_dk2 = dVariance_dk2(k1,k2,y(1:N-1),TimeStep);
dSigma_dk2 = diag(dSigma_dk2);

Sigma     = diag(Var);
InvSigma  = inv(Sigma);
dInvSigma_dk1 = -InvSigma*dSigma_dk1*InvSigma;
dInvSigma_dk2 = -InvSigma*dSigma_dk2*InvSigma;

L_deriv(1) = (dM_dk1'*InvSigma)*(y(2:end) - Mu) - 0.5*((y(2:end) - Mu)'*dInvSigma_dk1*(y(2:end) - Mu)) - 0.5*trace(InvSigma*dSigma_dk1);
L_deriv(2) = (dM_dk2'*InvSigma)*(y(2:end) - Mu) - 0.5*((y(2:end) - Mu)'*dInvSigma_dk2*(y(2:end) - Mu)) - 0.5*trace(InvSigma*dSigma_dk2);

G(1,1) = dM_dk1'*InvSigma*dM_dk1 + 0.5*trace(InvSigma*dSigma_dk1*InvSigma*dSigma_dk1);
G(1,2) = dM_dk1'*InvSigma*dM_dk2 + 0.5*trace(InvSigma*dSigma_dk1*InvSigma*dSigma_dk2);
G(2,1) = dM_dk2'*InvSigma*dM_dk1 + 0.5*trace(InvSigma*dSigma_dk2*InvSigma*dSigma_dk1);
G(2,2) = dM_dk2'*InvSigma*dM_dk2 + 0.5*trace(InvSigma*dSigma_dk2*InvSigma*dSigma_dk2);

J       = diag([(10^k1)*log(10) (10^k2)*log(10)]);
L_deriv = L_deriv*J;
G       = J'*G*J;

CurrentG = G;

% Inverse of G
CurrentInvG = inv(CurrentG);
    
CurrentFirstTerm = CurrentInvG*L_deriv';



for IterationNum = 1:NumOfIterations
        
    if mod(IterationNum,1000) == 0  && IterationNum < BurnIn
        disp([num2str(IterationNum) ' iterations completed.'])
        drawnow
    end
        
    %IterationNum
    ParasNew = Paras;
    Proposed = Proposed + 1;
    
    
    % Calculate the drift term
    Mean = ParasNew + (StepSize/(2))*CurrentFirstTerm';
    
    ParasNew = Mean + ( randn(1,D)*chol(StepSize*CurrentInvG) );
    
    % Calculate proposed Likelihood value
    
    Mu  = Calc_Mean(ParasNew(1),ParasNew(2),y(1:N-1),TimeStep);
    Var = Calc_Variance(ParasNew(1),ParasNew(2),y(1:N-1),TimeStep);
    %LogLikelihood  = sum( log( normpdf(y(2:end),Mu,sqrt(Var)) ) );
    LogLikelihood = pn(ParasNew(1), ParasNew(2), TimeStep, y);
    
    ProposedLJL    = LogLikelihood; % + LogPrior;


    ProbNewGivenOld = -sum(log(diag(chol(StepSize*CurrentInvG)))) - 0.5*(Mean-ParasNew)*(CurrentG/StepSize)*(Mean-ParasNew)';
    
    %%% Calculate probability of Old given New %%%
    
    % Calculate G based on new parameters
    dM_dk1     = dMean_dk1(ParasNew(1),ParasNew(2),y(1:N-1),TimeStep);
    dSigma_dk1 = dVariance_dk1(ParasNew(1),ParasNew(2),y(1:N-1),TimeStep);
    dSigma_dk1 = diag(dSigma_dk1);

    dM_dk2     = dMean_dk2(ParasNew(1),ParasNew(2),y(1:N-1),TimeStep);
    dSigma_dk2 = dVariance_dk2(ParasNew(1),ParasNew(2),y(1:N-1),TimeStep);
    dSigma_dk2 = diag(dSigma_dk2);

    Sigma     = diag(Var);
    InvSigma  = inv(Sigma);
    dInvSigma_dk1 = -InvSigma*dSigma_dk1*InvSigma;
    dInvSigma_dk2 = -InvSigma*dSigma_dk2*InvSigma;

    L_deriv(1) = (dM_dk1'*InvSigma)*(y(2:end) - Mu) - 0.5*((y(2:end) - Mu)'*dInvSigma_dk1*(y(2:end) - Mu)) - 0.5*trace(InvSigma*dSigma_dk1);
    L_deriv(2) = (dM_dk2'*InvSigma)*(y(2:end) - Mu) - 0.5*((y(2:end) - Mu)'*dInvSigma_dk2*(y(2:end) - Mu)) - 0.5*trace(InvSigma*dSigma_dk2);

    G(1,1) = dM_dk1'*InvSigma*dM_dk1 + 0.5*trace(InvSigma*dSigma_dk1*InvSigma*dSigma_dk1);
    G(1,2) = dM_dk1'*InvSigma*dM_dk2 + 0.5*trace(InvSigma*dSigma_dk1*InvSigma*dSigma_dk2);
    G(2,1) = dM_dk2'*InvSigma*dM_dk1 + 0.5*trace(InvSigma*dSigma_dk2*InvSigma*dSigma_dk1);
    G(2,2) = dM_dk2'*InvSigma*dM_dk2 + 0.5*trace(InvSigma*dSigma_dk2*InvSigma*dSigma_dk2);
    
    J       = diag([(10^ParasNew(1))*log(10) (10^ParasNew(2))*log(10)]);
    L_deriv = L_deriv*J;
    G       = J'*G*J;
    
    % Inverse of G
    InvG = inv(G);
    
    FirstTerm = InvG*L_deriv';
    
    % Calculate the drift term
    Mean = ParasNew + (StepSize/(2))*FirstTerm';
    
    try
        ProbOldGivenNew = -sum(log(diag(chol(StepSize*InvG)))) - 0.5*(Mean-Paras)*(G/StepSize)*(Mean-Paras)';
    catch
        %disp(ParasNew)
        ProbOldGivenNew = -1e300;
    end
    
    
    % Accept according to ratio
    Ratio = ProposedLJL + ProbOldGivenNew - CurrentLJL - ProbNewGivenOld;
        
        
    if Ratio > 0 || (Ratio > log(rand))
        CurrentLJL = ProposedLJL;
        
        CurrentG          = G;
        CurrentInvG       = InvG;
        CurrentFirstTerm  = FirstTerm;
        
        Paras = ParasNew;
        Accepted = Accepted + 1;
    end
        
    if mod(IterationNum, 100) == 1 && IterationNum < BurnIn
        Acceptance = Accepted/Proposed;
        disp(Acceptance)
        
        Proposed = 0;
        Accepted = 0;
    end
        
    % Save samples if required
    if IterationNum > BurnIn
        ParasSaved(IterationNum-BurnIn,:) = Paras;
        LJLSaved(IterationNum-BurnIn) = CurrentLJL;
    end
    
    % Start timer after burn-in
    if IterationNum == BurnIn
        disp('Burn-in complete, now drawing posterior samples.')
        tic;
    end
    
end

% Stop timer
TimeTaken = toc;

betaPosterior = ParasSaved;

CurTime = fix(clock);
save(['Results_mMALA_Simp_' num2str(floor(now)) '_' num2str(CurTime(4:6)) '.mat'], 'betaPosterior', 'LJLSaved', 'TimeTaken')


% Plot paths and histograms
figure(60)
plot(ParasSaved);
figure(61)
NumOfPlots = min(16,D);
for d = 1:NumOfPlots
    subplot(ceil(NumOfPlots/4),4,d)
    hist(ParasSaved(:,d))
end



end



function pn = pn(k1, k2, Time, n)

k1 = 10^k1;
k2 = 10^k2;

% Calculate the true probability

pn = 0;

for j = 1:length(n)-1
  x = n(j+1);
  n0 = n(j);
  t = Time(j);  
  k3 = k1/k2*(1-exp(-k2*t));
  value = 0;
  for i = 0:x 
      if n0 >= i
          value = value + nchoosek(n0, i)*exp(-k2*t*i)*(1-exp(-k2*t))^(n0-i)*poisspdf(x-i, k3);
      end
  end
  pn = pn + log(value);
  
end
  
end



% Moments
% k1 = immigration rate, k2 =death

function M = Calc_Mean(k1, k2, n0, t)

    k1 = 10^k1;
    k2 = 10^k2;

    % Moment Closure
    M = (k1-k1.*exp(-k2.*t)+n0.*exp(-k2.*t).*k2)./k2;
    
    % Crude
    %M = (k1-k2.*n0).*t;
end

function V = Calc_Variance(k1, k2, n0, t)

    k1 = 10^k1;
    k2 = 10^k2;
    
    % Moment Closure
    %V = (-n0.*exp(-2.*k2.*t).*k2+k1-k1.*exp(-k2.*t)+n0.*exp(-k2.*t).*k2)./k2;
    
    % Crude
    V = (k1+k2.*n0).*t;
end

% First-order 
function dM_dk1 = dMean_dk1(k1, k2,n0, t)

    k1 = 10^k1;
    k2 = 10^k2;
    
    % Moment Closure
    dM_dk1 = (1-exp(-k2.*t))./k2;
    
    % Crude
    %dM_dk1 = t;
end

function dM_dk2 = dMean_dk2(k1, k2,n0, t)

    k1 = 10^k1;
    k2 = 10^k2;
    
    % Moment Closure
    dM_dk2 = -k1.*(1-exp(-k2.*t))./k2.^2+k1.*t.*exp(-k2.*t)./k2-n0.*t.*exp(-k2.*t);
    
    % Crude
    %dM_dk2 = -n0.*t;
end

%Derivative of Variance wrt to theta
function dV_dk1 = dVariance_dk1(k1, k2,n0, t)

    k1 = 10^k1;
    k2 = 10^k2;
    
    % Moment Closure
    %dV_dk1 = (1-exp(-k2.*t))./k2;
    
    % Crude
    dV_dk1 = t;
end

function dV_dk2 = dVariance_dk2(k1, k2,n0, t)

    k1 = 10^k1;
    k2 = 10^k2;
    
    % Moment Closure    
    %dV_dk2 = (2.*n0.*t.*exp(-2.*k2.*t).*k2-n0.*exp(-2.*k2.*t)+k1.*t.*exp(-k2.*t)-n0.*t.*exp(-k2.*t).*k2+n0.*exp(-k2.*t))/k2-(-n0.*exp(-2.*k2.*t).*k2+k1-k1.*exp(-k2.*t)+n0.*exp(-k2.*t).*k2)./k2.^2;
    
    % Crude
    dV_dk2 = n0.*t;
end
