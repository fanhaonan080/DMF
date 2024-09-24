function[batch_splitted]=split2(batch,numsectors)

    noPointsWorker=diff(fix(linspace(0, size(batch,2), numsectors+1)));  % size of the splited alpha vectors: to be send to the workers
    batch_splitted=mat2cell(batch,size(batch,1),noPointsWorker);           % split alpha vector into  (~) even chunks




