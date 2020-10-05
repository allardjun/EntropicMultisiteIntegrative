% Integrative model for Lara

lambdaKArray = 2.^(0:0.25:1);
lambdaZArray = 2.^(-1:0.25:0);

paramArray = lambdaZArray;

kF0Array = logspace(-6,4,100);

PSSArray = zeros(numel(kF0Array),numel(paramArray));
ZSSArray = zeros(numel(kF0Array),numel(paramArray));

for iParam=1:numel(paramArray) % loop through parameter
    
    for iDose=1:numel(kF0Array) % loop through KP doses
        
        % biophysical parameters
        kK0 = 1; %s-1
        kF0 = kF0Array(iDose); %s-1
        kon0 = 0.1; %s-1
        koff = 0.01; %s-1
        
        N = 10; % number of sites
        
        lambdaK = 1;%lambdaKArray(iParam);%1;
        lambdaZ = lambdaZArray(iParam);%1;
        
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
    
    log(KPRatioArray);
    log(PSSArray);
    log(ZSSArray);
    
end % finished sweep through params

%%

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


