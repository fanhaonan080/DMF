function [delay_margin,DM] = DMF(A,B,tauinh,n)

g = 2e-3;    % grid step for theta
a = 2e-5;    % accuracy for the solution

theta_l=0; theta_u=2.1*pi;
points = ceil(abs(theta_l-theta_u)/g)+1;

p = gcp;                    % open new pool                    
poolsize = p.NumWorkers;

%% 
eg = eig(A+B);
phase = atan2(imag(eg),real(eg));
phase = phase + (phase<0)*2*pi;
om = abs(eg);
id_nonzero = om>1e-10;
tauinhDM = min((phase(id_nonzero)-pi/2)./om(id_nonzero));
% tauinhDM = pi/2;
if any((abs(phase(id_nonzero))-pi/2)<0)
    error('Delay-free system is unstable')
end
if tauinh > tauinhDM
    error('Inherent delay (tauinh) is too large')
end

%%
%==========================================================================
% Split theta-grid (coarse) in chunks
%==========================================================================
theta = linspace(theta_l,theta_u,points);
[theta_splitted]=split2(theta,poolsize);

delta = abs(theta_l-theta_u)/(points-1);

theta_chunks=theta_splitted;

for i=1:poolsize-1
    theta_chunks{i}=[theta_splitted{i},theta_splitted{i}(end)+delta];   % extended intervals 
end


%%
lambdaInter{poolsize,1}=[];
thetaInter{poolsize,1}=[];

parfor ii=1:poolsize
[FamLambda,theta] = numbering2(n,theta_chunks{ii},A,B);
ReFamLambda = real(FamLambda);
ImFamLambda = imag(FamLambda);

omega = sqrt(ReFamLambda.^2+ImFamLambda.^2);

obj = ReFamLambda.*cos(omega*tauinh) + ImFamLambda.*sin(omega*tauinh);


% theta-solutions in the coarse grid
%----------------------------------------------------------------------
indices_zx = zeros(n,n);
indices = find(sum(abs(obj'))>g); % rule out lambda=0 eigenval

for kk = 1:length(indices)
    zx = zci(obj(indices(kk),:));
    indices_zx(1:length(zx),indices(kk)) = zx;
end

[q,m,index]=find(indices_zx);   % q-th theta solution for the coarse grid
% m-th factor
if isempty(q)
    MinLocalTau(ii)=NaN; % jump to the "extraction of DM" part of the code
    MinLocalTheta(ii)=NaN;
else
    thetaInter{ii} = [theta(index)',theta(index+1)'];                % theta interval were the candidate solution may be found
    lambdaInter{ii} = [diag(FamLambda(m,index)),diag(FamLambda(m,index+1))];  % eigenvalue interval were the candidate solution may be found
end
end

%%
%==========================================================================
% Split theta(interval)-solutions and the corresponding lambda values
%==========================================================================
thetaInter=thetaInter(~cellfun('isempty',thetaInter));   %rule out empty cells
lambdaInter=lambdaInter(~cellfun('isempty',lambdaInter));

thetaInter=cell2mat(thetaInter).';
lambdaInter=cell2mat(lambdaInter).';

nn=min(poolsize,size(thetaInter,2)); % with q empty, then nn=0


[thetaInter_splitted]=split2(thetaInter,nn);
[lambdaInter_splitted]=split2(lambdaInter,nn);

%==========================================================================
% Second batch in parallel
%==========================================================================
parfor ii=1:nn
    % increase accuracy of theta-solutions
    %----------------------------------------------------------------------
    thetaCRIT=zeros(size(thetaInter_splitted{ii},2),1);
    lambdaCRIT=zeros(size(thetaInter_splitted{ii},2),1);
    for jj=1:size(thetaInter_splitted{ii},2)
        [thetaCRIT(jj,:),lambdaCRIT(jj,:)] = intersections2(a,thetaInter_splitted{ii}(1,jj), thetaInter_splitted{ii}(2,jj),lambdaInter_splitted{ii}(1,jj),lambdaInter_splitted{ii}(2,jj),A,B,tauinh);
    end
    
    % (omega,tau)-solutions in the refined grid
    %----------------------------------------------------------------------
    ReLambdaCRIT = real(lambdaCRIT);
    ImLambdaCRIT = imag(lambdaCRIT);
    Omega = sqrt(ReLambdaCRIT.^2+ImLambdaCRIT.^2);
    Theta = thetaCRIT;
    
%     indexOm = find(Omega); % consider only non-zero omegas
%     Omega = Omega(indexOm);
%     Theta = Theta(indexOm);
    
    indexTh = find(Theta>1e-04);
    Omega = Omega(indexTh);
    Theta = Theta(indexTh);
    
    %     Numerator = atan2(sqrt(1-Theta.^2),Theta);
    %     indexNum = find(Numerator<0);
    %     Numerator(indexNum) = Numerator(indexNum)+2*pi;
    %     Tau = -Numerator./Omega;
    Tau = Theta./Omega;
    %     indexTau = find(Tau<0);
    %     Tau(indexTau) = Tau(indexTau)+2*pi./Omega(indexTau);
    Tau(Tau==0) = NaN; % rule out tau=0
    
      
    [MinTau,indexMinTau]=min(Tau);
    MinTheta=Theta(indexMinTau);
    
%     MinLocalTau(ii)=MinTau;
%     MinLocalTheta(ii)=MinTheta;
    if isempty(indexMinTau)
        MinLocalTau(ii)=NaN;
        MinLocalTheta(ii)=NaN;
    else
        MinLocalTau(ii)=MinTau;
        MinLocalTheta(ii)=MinTheta;
    end
end

%==========================================================================
% extraction of DM
%==========================================================================
[DM,indexDM]=min(MinLocalTau(MinLocalTau>0),[],'omitnan'); % minimum excluding NaN and zero values

% ThetaDM=MinLocalTheta(indexDM);

indexDM(isnan(DM))=[];                     % in case MinLocalTau contains only NaN values then make indexDm empty
if isempty(indexDM)
    delay_margin = 'infinite';
    return
end

[delay_margin] = num2str(DM);

return


