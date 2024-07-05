function [data_n, y, u] = normalize_data(data,range)
%  It is necessary to normalize. Get the normalization according to the
%  complete dataset
[~,y.ctr,y.scl] = normalize(cell2mat({data.y}'),"range",[-1,1]);
% the same for the input
if isfield(data,"u")
  [~,u.ctr,u.scl] = normalize(cell2mat({data.u}'),"range",range);
end
% apply the normalization to every sample. according t the center and the
% scale
if isfield(data,'u')
  data_n = arrayfun(@(x) struct('y', normalize(x.y,"center",y.ctr, ...
  	"scale",y.scl), ...
  	'u', normalize(x.u,"center",u.ctr, ...
  	"scale",u.scl), ...
  	't', x.t),data); % keep the time as is
else
  data_n = arrayfun(@(x) struct('y', normalize(x.y,"center",y.ctr,...
    "scale",y.scl),...
    't',x.t,data));
end

end
