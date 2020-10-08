% Integrative model for Lara

lambdaKArray = 2.^(0:0.1:1);
lambdaZArray = 2.^(-2:0.1:0);

paramArray = lambdaZArray;

kF0Array = logspace(-5,5,400);

PSSArray = zeros(numel(kF0Array),numel(paramArray));
ZSSArray = zeros(numel(kF0Array),numel(paramArray));


nH_Z_Array = zeros(0,numel(paramArray));
nH_P_Array = zeros(0,numel(paramArray));

figure(4);clf(4);
figure(5);clf(5);

for iParam=1:numel(paramArray) % loop through parameter
    
    for iDose=1:numel(kF0Array) % loop through KP doses
        
        % biophysical parameters
        kK0 = 1; %s-1
        kF0 = kF0Array(iDose); %s-1
        kon0 = 0.1; %s-1
        koff = 0.01; %s-1
        
        N = 10; % number of sites
        
        lambdaK = 1.3;%lambdaKArray(iParam);%1;
        lambdaZ = lambdaZArray(iParam);%1;
        
        % system
        
        dnPdt =@(nP,nZ) kK0*lambdaK.^nP.*lambdaZ.^nZ.*(N-nP) - kF0*lambdaK.^nP.*lambdaZ.^nZ.*(nP-nZ);
        dnZdt =@(nP,nZ) kon0*lambdaK.^nP.*lambdaZ.^nZ.*(nP-nZ) - koff*nZ;
        
        %dnPdt =@(nP,nZ) kK0*lambdaK.^nP.*lambdaZ.^nZ.*(N-nP) - kF0*(nP-nZ); % CONSTANT DEPHOSPHORYLATION RATE
        
        [T,X] = ode23s( @(t,x)[dnPdt(x(1),x(2));dnZdt(x(1),x(2))], [0,1e3], [0,0]);
        
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
                
        [nHill_finiteDifference, indexMax] = max(diff(log(PSSArray(:,iParam)./(N-PSSArray(:,iParam))))./diff(log(KPRatioArray')));
        
        if nHill_finiteDifference > 1.001
            domainStart = max(1,indexMax-10);
            if indexMax == 1
                domainEnd = 6;
            else
                domainEnd = indexMax+10;
            end
            
            y = log(PSSArray(:,iParam)./(N-PSSArray(:,iParam)));
            x = log(KPRatioArray');
            diffx = diff(x(domainStart:domainEnd));
            diffy = diff(y(domainStart:domainEnd));
            slope = diffy./diffx;
            
            fit = polyfit(x(domainStart:domainEnd-1),slope,3);
            slope_fit = polyval(fit,x(domainStart:domainEnd-1));
            HillCoeffMaxSlope = max(slope_fit);
            
            nH_P_Array(iParam) = HillCoeffMaxSlope;
            
            figure(5); hold on;
            subplot(4,1,1);hold on;
            plot(KPRatioArray,PSSArray(:,iParam))
            set(gca,'xscale','log');
            
            subplot(4,1,2);hold on;
            plot(KPRatioArray,PSSArray(:,iParam)./(N-PSSArray(:,iParam)))
            set(gca,'xscale','log','yscale','log');
            
            subplot(4,1,3);hold on;
            plot(x,y)
            
            subplot(4,1,4);hold on;
            plot(x(1:end-1),diff(y)./diff(x))
         
            
            figure(4); hold on;
            plot(x(domainStart:domainEnd-1), slope,'s')
            plot(x(domainStart:domainEnd-1), slope_fit,'-');
            
            drawnow;
            
        else
            nH_P_Array(iParam) = 1;
        end
               
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
%set(gca,'yscale','log');
ylabel('Phosphorylation Hill coefficient');
xlabel('lambda')

subplot(2,1,2); hold on; box on;
plot(paramArray, nH_Z_Array,'d-');
set(gca,'xscale','log');
%set(gca,'yscale','log');
ylabel('ZAP70 binding Hill coefficient');
xlabel('lambda');

%set(gca,'ylim', [-1,10])
