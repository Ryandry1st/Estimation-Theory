% Multiple-model estimation example: The solar panel deployment problem
% User defined
clear all;
close all;

sims = 400;
xtilde = zeros(sims, 1000, 2);
plotloc = [400, 400, 1800, 1300];

mincost = 47.4;
ideal = 0.42378;
format long;

% load provided data
load('mmExampleModeSwitchingData.mat')

% uncomment to test consistency over a range
% for lowerbound=0.4237:0.00002:0.439
    zero_mean = 0;
    lower_bound = ideal;
    switching = 2;
    switchtime = 20;
    % set to 1 to include transistion probability
    transition = 0;

    %----- Simulation parameters
    dt = 0.1;
    Nsim = 1000;
    % tkhist = [0:Nsim-1]'*dt;
    % nmod is the number of models considered
    nmod = 3;  
    nx = 2;
    nz = 1;
    nu = 1;

    %----- Random number seed
    % n_seed = round(sum(clock*100));
    % For reproducible results
    n_seed = 1;
    randn('state', n_seed);


    for expectationrun=1:sims
    %----- Storage matrices
    % xtruekhist = zeros(Nsim,nx);
    % ztruekhist = zeros(Nsim,nz);
    xhatkhist = zeros(Nsim,nx);
    Pkhist = zeros(Nsim,nx);    
    % xhatkMhist is an array of all the mode-conditioned estimates
    xhatkMhist = zeros(Nsim,nx,nmod); 
    PkMhist = zeros(Nsim,nx,nx,nmod);
    mukhist = zeros(Nsim,nmod);
    nukMhist = zeros(Nsim,nz,nmod);
    LambdakMhist = zeros(Nsim,nmod);
    FkM = zeros(nx,nx,nmod);
    GkM = zeros(nx,nu,nmod);
    invSkMhist = zeros(Nsim, nmod);

    %----- Set up multiple system models
    cSet = [0.1;0.5;1];  hzSet = [1;2;3];
    avec = cSet./hzSet;
    bvec = 1./hzSet;
    for j=1:nmod
      a = avec(j);
      b = bvec(j);
      Amat = [0, 1; 0, -a];
      % Discrete-time state transition matrix
      % FkM(:,:,j) = Amat*dt+eye(nx, nx);
      FkM(:, :, j) = expm(Amat*dt);
      % Discrete-time control matrix
      GkM(:,:,j) = [b/a^2*(exp(-a*dt) + a*dt - 1); -(b/a*(exp(-a*dt) - 1))];
    end

    %----- Truth model
    % switchHist gives the model switching time history and switchIndex gives
    % the indices at and after which each model obtains.  For the static case,
    % switchHist is either 1, 2, or 3, and switchIndex = 1;
    if switching == 1
        switchHist = randi([1, 3], 1, Nsim*dt/switchtime); 
        switchIndex = 1:switchtime/dt:Nsim;
        % Normally I would randomly assign states
        % but since it impacts the dynamics equations,
        % a good estimate cannot be had without some knowlege of the dynamics
        % at that time. So I attempt to determine them and encode them.
    elseif switching == 2
        switchHist = [1, 3, 1];
        switchIndex = [1, 200, 800];
    else
        switchHist = [1];
        switchIndex = [1];
    end
    Qk = 0.001*diag([0.1 1]); Rq = chol(Qk);
    vtruekhist = (Rq'*randn(nx,Nsim))';
    Rkp1 = 0.1; Rr = chol(Rkp1);
    wtruekhist = (Rr'*randn(nz,Nsim))';
    % utruekhist = 2*randn(Nsim,nu);
    Hkp1 = [1 0];
    x1 = [0;0.1];

    %----- Generate truth-model states and measurements
    xk = x1;
    ii = 1;
    for k=1:Nsim-1
      if k==switchIndex(ii)
        Fk = FkM(:,:,switchHist(ii));
        Gk = GkM(:,:,switchHist(ii));
        if ii<length(switchHist)
          ii = ii+1;
        end
      end
      kp1 = k + 1;
      uk = utruekhist(k,:)';
      vk = vtruekhist(k,:)';
      wkp1 = wtruekhist(kp1,:)';
      xkp1 = Fk*xk + Gk*uk + vk;
      zkp1 = Hkp1*xkp1 + wkp1;
      % xtruekhist(k,:) = xk';
      % ztruekhist(kp1,:) = zkp1';
      xk = xkp1;
    end

    %----- Run the Multiple-model filter
    % Set up initial states and covariances
    P1 = 10*eye(nx);
    xhat1 = x1 + (chol(P1))'*randn(nx,1);
    for j=1:nmod
      PkMhist(1,:,:,j) = P1;
      % Initialize the mode probabilities as equally probable
      mukhist(1,j) = 1/nmod;
      xhatkMhist(1,:,j) = xhat1';
    end
    for k=1:Nsim-1
      kp1 = k + 1;
      zkp1 = ztruekhist(kp1,:)'; 
      uk = utruekhist(k,:)';
      for j=1:nmod
        % Propagation step
        Fkj = FkM(:,:,j);
        Gkj = GkM(:,:,j);
        xhatkj = xhatkMhist(k,:,j)';
        Pkj = squeeze(PkMhist(k,:,:,j));
        xbarkp1j = Fkj*xhatkj + Gkj*uk;
        Pbarkp1j = Fkj*Pkj*Fkj' + Qk;

        % Measurement update
        nukp1j = zkp1 - Hkp1*xbarkp1j;         
        Skp1j = Hkp1*Pbarkp1j*Hkp1' + Rkp1;
        invSkp1j = inv(Skp1j);                 
        Wkp1j = Pbarkp1j*Hkp1'*invSkp1j;
        xhatkp1j = xbarkp1j + Wkp1j*nukp1j;
        Pkp1j = Pbarkp1j - Wkp1j*Skp1j*Wkp1j';

        % Calculate the likelihood of the current innovation
        Lambdakp1j = normpdf(nukp1j, zero_mean, Skp1j);

        % Store state estimate and covariance, innovation, and likelihood
        xhatkMhist(kp1,:,j) = xhatkp1j';
        PkMhist(kp1,:,:,j) = Pkp1j;
        nukMhist(kp1,:,j) = nukp1j';
        LambdakMhist(kp1,j) = Lambdakp1j;
        % also store the inverse S for consistency testing
        invSkMhist(kp1, j) = invSkp1j;
      end


      % Update the mode probabilities
      mukvec = mukhist(k,:)';
      Lambdakp1vec = LambdakMhist(kp1,:)';
      % mukp1vec is the nmod-by-1 vector that contains the mode probabilities
      % for the nmod modes at time index kp1; it must satisfy sum(mukp1vec) = 1.
      mukp1vec = Lambdakp1vec.*mukvec./sum(Lambdakp1vec.*mukvec);

      % Ad hoc adjustment
      % assume low probability of switching given previous state due to
      % frequentist approach
      % most likely previous mode
      if size(switchHist, 2) > 1
          for anothervariable=1:nmod
              if mukp1vec(anothervariable) < lower_bound
                  mukp1vec(anothervariable) = lower_bound;
              end
          end
          if transition == 1
              [largestmu, modesel] = max(mukhist(k, :));
              pij = [1; 1; 1]/40;
              pij(modesel) = 19/20;
              mukvec = mukvec.*pij;
              mukp1vec = Lambdakp1vec.*mukvec./sum(Lambdakp1vec.*mukvec);
          end
          % renormalize
          mukp1vec = mukp1vec./(sum(mukp1vec));
      end

      mukhist(kp1,:) = mukp1vec';

      % Calculate the combined state estimate and covariance
      xMdum = reshape(squeeze(xhatkMhist(kp1,:,:)),nx,nmod);
      mudum = mukhist(kp1,:)';
      % Calculate the MMSE estimate averaged over all the modes
      xhatkp1 = (xMdum*mudum);
      % Calculate the estimation error covariance of xhatkp1
      Pkp1 = zeros(nx,nx);
      for j=1:nmod
        Pkp1j = squeeze(PkMhist(kp1,:,:,j));
        xhatkp1j = xhatkMhist(kp1,:,j)';
        mukp1j = mukhist(kp1,j);
        Pkp1 = Pkp1 + mukp1j*(Pkp1j + (xhatkp1j-xhatkp1)*(xhatkp1j-xhatkp1)');
      end
      xhatkhist(kp1,:) = xhatkp1'; 
      Pkhist(kp1,:) = diag(Pkp1)';  
    end

    xtilde(expectationrun, :, :) = xtruekhist - xhatkhist;
    for thisstep=2:Nsim
        ev(expectationrun, thisstep) = squeeze(xtilde(expectationrun, thisstep, :))'*inv(diag(squeeze(Pkhist(thisstep, :))))*squeeze(xtilde(expectationrun, thisstep, :));
    end
    % record the expectation for each time step, as more simulations are added
    % to the averaging
    meanev(expectationrun, :) = mean(ev(1:expectationrun, :), 1);
    end

    % determine how consistent the filter was
    cost = abs(sum(meanev(end, :)-nx))
    if cost < mincost
        mincost = cost
        ideal = lower_bound
        bestev(:, :) = meanev(:, :);
    end

% end
%% 
%----- Display results
figure(1);clf;
plot(tkhist,mukhist);shg;
xlabel('Time (s)');
ylabel('Probability');
legend('\mu_1', '\mu_2','\mu_3');
title('Mode probability time histories');
set(gcf, 'position', plotloc);
ylim([0.1 0.8]);
xlim([0, max(tkhist)]);

figure(2);clf;
subplot(211)
iidum = 1:Nsim-1;
plot(tkhist(iidum), xtruekhist(iidum,1) - xhatkhist(iidum,1));
hold on;
plot(tkhist(iidum), sqrt(Pkhist(iidum,1)), 'k');
plot(tkhist(iidum), -sqrt(Pkhist(iidum,1)), 'k');
legend('Theta difference', 'Error Cov');
ylabel('\Delta \theta_z (rad)');
title('Estimation errors and covariances');
xlim([0, max(tkhist)]);
set(gcf, 'position', plotloc);


subplot(212)
iidum = 1:Nsim-1;
plot(tkhist(iidum), xtruekhist(iidum,2) - xhatkhist(iidum,2));
hold on;
plot(tkhist(iidum), sqrt(Pkhist(iidum,2)), 'k');
plot(tkhist(iidum), -sqrt(Pkhist(iidum,2)), 'k');
legend('Omega difference', 'Error Cov');
ylabel('\Delta \omega_z (rad/s)');
xlabel('Time (s)');
xlim([0, max(tkhist)]);
set(gcf, 'position', plotloc);
