%
% [face_labels] = get_face_labels()
%
% output:
%   A 1xN cell array, where N is the number of features.

function [face_labels] = get_face_labels()
    
    face_labels = {};
    types = {'', 'delta: ', 'double_delta: '};
    
    for i=1:length(types)
        elem = types{i};
        face_labels = [face_labels [elem 'mean'] [elem 'abs mean']];
        face_labels = [face_labels [elem 'max'] [elem 'abs max'] ...
        [elem, 'min'] [elem 'abs min']];
        face_labels = [face_labels [elem 'max+min'] [elem 'max-min']];
        face_labels = [face_labels [elem 'sixth power mean']];
        face_labels = [face_labels [elem 'std'] [elem 'var']];
        face_labels = [face_labels [elem 'kurtosis'] [elem 'skewness']];
        face_labels = [face_labels [elem 'sum'] [elem 'median']];
    end
end