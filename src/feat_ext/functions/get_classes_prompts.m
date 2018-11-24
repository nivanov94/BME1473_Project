%
% [classes] = get_classes_prompts(prompts)
%
% input:
%   prompts     - the prompts struct.
% output:
%   the classes, yo.
function [class_labels] = get_classes_prompts(prompts)
    
    classes = struct();
    classes.iy  =   0;
    classes.uw  =   1;
    classes.m   =   2;
    classes.n   =   3;
    classes.piy =   4;
    classes.tiy =   5;
    classes.diy =   6;
    classes.gnaw =  7;
    classes.knew =  8;
    classes.pat =   9;
    classes.pot =   10;
    
    class_labels = [];

    for i=1:length(prompts)
      prompt = regexprep(prompts{i}, '/', '');
      class_labels = [class_labels; classes.(prompt)];
    end
    
end