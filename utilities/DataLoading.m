function s = DataLoading(graph_select,normalizeAB,Randomize,savdir,n)
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
% 
% if graph_select == 1
%     if normalizeAB == 0
%         filename = sprintf('Fully_tau_scan_n%d.mat',n);
%     else
%         filename = sprintf('Normalized_Fully_tau_scan_n%d.mat',n);
%     end
% elseif graph_select == 2
%     if normalizeAB == 0
%         filename = sprintf('BilateralRing_tau_scan_n%d.mat',n);
%     else
%         filename = sprintf('Normalized_BilateralRing_tau_scan_n%d.mat',n);
%     end
% elseif graph_select == 3
%     if normalizeAB == 0
%         filename = sprintf('DirectedRing_tau_scan_n%d.mat',n);
%     else
%         filename = sprintf('Normalized_DirectedRing_tau_scan_n%d.mat',n);
%     end
% end

s = load(fullfile(savdir,filename));

end