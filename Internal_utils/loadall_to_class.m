function neuron = loadall_to_class(neuron, savedfile)

    if isstruct(savedfile)
        tempvar = savedfile; 
    else 
        [~,~,fmt] = fileparts(savedfile);
        if strcmpi(fmt, '.mat')
            tempvar = load(savedfile);
        else
            error('Unrecognized file format'); 
        end 
    end
    
    props = properties(neuron);  
    for p=1:numel(props)
        neuron.(props{p}) = tempvar.neuron.(props{p}); 
    end
    
    flds = fieldnames(tempvar);
    m=1;
    while m<=numel(flds) 
        if isempty(intersect(flds{m},props)) && ~strcmpi(flds{m},'neuron')
            assignin('base',flds{m},tempvar.(flds{m}))   
        end
        m = m+1; 
    end
    
end