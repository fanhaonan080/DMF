close all
clear
addpath('.\utilities')
if 0==exist('ScanningResult','dir')       % check if Data folder already exists
    mkdir ScanningResult                  % create new directory to save the data
end
savdir = '.\ScanningResult';
%% 
n_list = [4];%[5 6 10 50 100 200];        % system size
graph_select = 3;   % 1 strongly connected, 2 bilateral ring, 3 directed ring
normalizeAB = 1;
iter = 100;        % number of tau_inh
rand_seed = 10; %randi(100,1);
Randomize_list = [0];    % 0:original matrices, 1:plus 0~10percent, 2:minus 0~10percent
%% 
Plot_Result = 1;
savedata = 1;
%%
for ii_rand = 1:length(Randomize_list)
    Randomize = Randomize_list(ii_rand);
for i = 1:length(n_list)
    n = n_list(i)
    [A,B] = MatrixAB_Generation(graph_select,normalizeAB,n);
%     [A,B] = matrices(n);
    %% Randomize A and B 
    rng(rand_seed)
    seed = randi(100,1,2);
    if Randomize == 1
        % plus 10 percent
        rng(seed(1))
        RM_plus = rand(n,n);   % random matrix
        scaling_factor_plus10 = 1+RM_plus*0.1;
        B = B.*scaling_factor_plus10;
        A = -diag(sum(B,2));
    elseif Randomize == 2
        % minus 10 percent
        rng(seed(2))
        RM_minus = rand(n,n);   % random matrix
        scaling_factor_minus10 = 1-RM_minus*0.1;
        B = B.*scaling_factor_minus10;
        A = -diag(sum(B,2));
    elseif Randomize == 0
    else
        error('Randomize should be 0, 1 or 2.')
    end
    
    
    %% tau_inh Margin
    tauMargin = tau_inh_finding(A+B);
    %% DM-finding
    p = -min(eig(A));
    tau_min = 1/2/p;
%     tauinh = linspace(0.5,tauMargin,iter);
%     tauinh = linspace(0,tauMargin,iter);
    tauinh = linspace(tau_min,tauMargin,iter);
    store = zeros(length(tauinh),2);
    cal_time = zeros(length(tauinh),1);
%     tic
    for ii = 1:length(tauinh)
        tic
        [delay_margin,DM] = DMF(A,B,tauinh(ii),n);%DMF(tauinh(ii),n);
        if isempty(DM)
            DM = Inf;
        end
        store(ii,:) = [tauinh(ii),DM];
        cal_time(ii,:) = toc;
    end
%     cal_time = toc
    
    %% Data storing
    if savedata == 1
        filename = DataSaving(rand_seed,store,n,cal_time,graph_select,normalizeAB,Randomize,savdir);
    end
end
end
%% Plot
if Plot_Result == 1
    figure,clf,hold on
    for jj = 1:length(n_list)
        s = DataLoading(graph_select,normalizeAB,Randomize,savdir,n_list(jj)); % load 'store' and 'n'
        line_name = sprintf('n = %d',n_list(jj));
        plot(s.store(:,1),s.store(:,2),'x','DisplayName',line_name)
    end
    grid
    legend
    xlim([0,max(s.store(:,1))*1.1])
    % ylim([-0.05,10])
    ylabel('\tau_0^\ast'), xlabel('\tau_{inh}')
    title('\tau_0^\ast V.S. \tau_{inh}')
end

%%
% system('shutdown -s')
