function [feats, p] = featSel_corr( in, class, N )
%  featSel_corr 
%
%   'in' is NxD
%       feats is a vector of the  N 'top' indices
  assert( length(class) == size(in, 1), 'number of observations and classes differ' )

  class = reshape(class,length(class),1);
%  class = repmat(class, size(in,2), 1);
%  class = reshape(class, size(in,2)*size(in,3), 1);
%  in = reshape(in, [size(in,2)*size(in,3), size(in,1)]);

  nonzerocols = find(any(in));
  b = zeros( 1, length(class));
  for i=nonzerocols%1:size(in,2)
    bb = corrcoef(in(:,i), double(class));
    b(i) = bb(2);
  end
  %b = corrcoef( [in double(class)]); b = b(1:(end-1),end);
  [p,i] = sort(abs(b), 'descend');
  feats = i(1:min(N,length(i)));
  p     = p(1:min(N,length(p)));

  return