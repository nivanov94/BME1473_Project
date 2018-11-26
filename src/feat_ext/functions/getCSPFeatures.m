function features = getCSPFeatures(E,W,band,Fs)

[z,p,k] = butter(10,band/(Fs/2));
sos = zp2sos(z,p,k);

E_filt = sosfilt(sos,E,2); % filter the readings from each electrode from each sample

Z = W*E_filt;

% now calculate the variance of each row and generate scores
v = var(Z,0,2);
features = log10(v/sum(v));
