% Integrative model for Lara

lambdaKArray = 2.^(0:0.1:1);
lambdaZArray = 2.^(-1:0.1:0);

paramArray = lambdaKArray;

kF0Array = logspace(-6,4,100);

PSSArray = zeros(numel(kF0Array),numel(paramArray));
ZSSArray = zeros(numel(kF0Array),numel(paramArray));


nH_Z_Array = zeros(0,numel(paramArray));
nH_P_Array = zeros(0,numel(paramArray));


for iParam=1:numel(paramArray) % loop through parameter
    
    for iDose=1:numel(kF0Array) % loop through KP doses
        
        % biophysical parameters
        kK0 = 1; %s-1
        kF0 = kF0Array(iDose); %s-1
        kon0 = 0.1; %s-1
        koff = 0.01; %s-1
        
        N = 10; % number of sites
        
        lambdaK = lambdaKArray(iParam);%1;
        lambdaZ = 1;%1;%lambdaZArray(iParam);%1;
        
        % system
        
        dnPdt =@(nP,nZ) kK0*lambdaK.^nP.*lambdaZ.^nZ.*(N-nP) - kF0*lambdaK.^nP.*lambdaZ.^nZ.*(nP-nZ);
        dnZdt =@(nP,nZ) kon0*lambdaK.^nP.*lambdaZ.^nZ.*(nP-nZ) - koff*nZ;
        
        %dnPdt =@(nP,nZ) kK0*lambdaK.^nP.*lambdaZ.^nZ.*(N-nP) - kF0*(nP-nZ); % CONSTANT DEPHOSPHORYLATION RATE
        
        [T,X] = ode15s( @(t,x)[dnPdt(x(1),x(2));dnZdt(x(1),x(2))], [0,1e3], [0,0]);
        
        if(0)
            figure(1); %clf;
            hold on;
            plot(T,X(:,1),'-', 'color', [0.5 0 1]); % purple for phosphorylated
            plot(T,X(:,2),'-r'); % red for ZAP
            
            display(iDose);
        end
        
        PSSArray(iDose,iParam) = X(end,1);
        ZSSArray(iDose,iParam) = X(end,2);
        
    end % finished loop through KP doses
    
    % Compute Hill
    KPRatioArray = kK0./kF0Array;
    
    if(any(diff(PSSArray(:,iParam))>0.001))
        nH_P_Array(iParam) = -1; % set to -1 if it's nonmonotonic increasing.
    else
        nH_P_Array(iParam) = max(diff(log(PSSArray(:,iParam)./(N-PSSArray(:,iParam))))./...
            diff(log(KPRatioArray')));
            
            
% legacy code from polymer-c analysis
%             % compute slope
%             diffy = diff(log10((avgSteadyState(domainStart:domainEnd))./(1-avgSteadyState(domainStart:domainEnd))));
%             diffx = diff(log10(kinaseIntrinsicRate(domainStart:domainEnd)));
%             slope = diffy./diffx;

%             % previous, alternate methods
%             %diffy = diff(movmean(log10((avgSteadyState(domainStart:domainEnd))./(1-avgSteadyState(domainStart:domainEnd))),3));
%             %diffy = diff((log10((avgSteadyState(domainStart:domainEnd))./(1-avgSteadyState(domainStart:domainEnd)))));

%             % fit slopes to cubic
%             fit = polyfit(log10(kinaseIntrinsicRate(domainStart+1:domainEnd)),slope,3);
%             slope_fit = polyval(fit,log10(kinaseIntrinsicRate(domainStart+1:domainEnd)));

%             % calculate sum of squared residuals and rmse
%             SSR = sum((slope_fit-slope).^2)
%             slope_rmse = sqrt(SSR./length(slope))

%             % set max slope and error
%             HillCoeffMaxSlope = max(slope_fit);
%             HillCoeffMaxSlope_std = slope_rmse;


    end
    
    if(any(diff(ZSSArray(:,iParam))>0.001))
        nH_Z_Array(iParam) = -1;  % set to -1 if it's nonmonotonic increasing.
    else
        nH_Z_Array(iParam) = max(diff(log(ZSSArray(:,iParam)./(N-ZSSArray(:,iParam))))./...
            diff(log(KPRatioArray')));
    end
    
    
    
end % finished sweep through params

%%

% Dose response curves

figure(2); clf;
subplot(2,1,1); hold on; box on;
plot(KPRatioArray, PSSArray);
set(gca,'xscale','log');
xlabel('Kinase Phosphatase ratio (ratio of intrinsic rate)');
ylabel('Number phosphorylation (out of 10)')

legend(num2str(paramArray','%3.2f'),'location','southeast')

subplot(2,1,2); hold on; box on;
plot(KPRatioArray, ZSSArray);
set(gca,'xscale','log');
ylabel('Number of ZAP70 bound (out of 10)')
xlabel('Kinase Phosphatase ratio (ratio of intrinsic rate)');

% Hill coefficients

figure(3); %clf;
subplot(2,1,1); hold on; box on;
plot(paramArray, nH_P_Array,'d-');
set(gca,'xscale','log');
ylabel('Phosphorylation Hill coefficient');
xlabel('lambda')

subplot(2,1,2); hold on; box on;
plot(paramArray, nH_Z_Array,'d-');
set(gca,'xscale','log');
ylabel('ZAP70 binding Hill coefficient');
xlabel('lambda');

set(gca,'ylim', [-1,10])
