% [prompts] = get_prompts(labels_fn);
%
% input:
%   labels_fn       - The filename of the file containing all the prompts.
% output:
%   prompts         - A 1xN cell array, where element i is the name of prompt i.
function [prompts] = get_prompts(labels_fn)
    
    prompts = {};
    fid = fopen(labels_fn);
    p = fgets(fid);
    
    while p ~= -1
        p = strtrim(p);
        prompts = [prompts p];
        
        p = fgets(fid);
    end
    fclose(fid);
end