function [C1, C2] = getCSPCompCovMat(E,labels,normalized_passband,ICA_Xform)

% begin by filter the each epoch
[z,p,k] = butter(10,normalized_passband);
sos = zp2sos(z,p,k);

Efilt = cell(1,length(E));
for e = 1:length(E)
  Efilt{e} = sosfilt(sos,ICA_Xform*double(E{e}),2);
end

[C1, C2] = getAverageClassCovMats(Efilt,labels);
