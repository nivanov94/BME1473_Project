function norm = normalize_feats( data ) 
%Vectorize later
    d_offsets = zeros(1,size(data,2));
    d_extents = zeros(1,size(data,2));
    norm = zeros( size(data) );
    for d=1:size(data,2)
      d_offsets(d)= min(data(:,d));
      norm(:,d) = data(:,d)-d_offsets(d);
      d_extents(d)= max(norm(:,d));
      norm(:,d) = norm(:,d)./d_extents(d);
    end
return