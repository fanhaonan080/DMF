function [filename] = DataSaving(rand_seed,store,n,cal_time,graph_select,normalizeAB,Randomize,savdir)

type = sprintf('tau_scan_n%d',n);
if graph_select == 1
    graph = 'Fully_'; 
elseif graph_select == 2
    graph = 'BilateralRing_';
elseif graph_select == 3
    graph = 'DirectedRing_'; 
end

if normalizeAB == 0
    normal = '';
else
    normal = 'Normalized_';
end

if Randomize == 0
    final = '.mat';
elseif Randomize == 1
    final = '_plus.mat';
elseif Randomize == 2
    final = '_minus.mat';
else
    error('Randomize should be 0, 1 or 2.')
end  

filename = strcat(normal,graph,type,final);

save(fullfile(savdir,filename),'store','n','cal_time','rand_seed');

end