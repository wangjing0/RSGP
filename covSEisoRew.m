function  [varargout] =covSEisoRew(varargin)
% covariance function: squared exponential, isotropic, and reward dependent
% Jing Wang   jingwang.physics(a)gmail.com
varargout = cell(max(1,nargout),1);
[varargout{:}] = covScaleRew({'covSE','iso',[]},varargin{:});

