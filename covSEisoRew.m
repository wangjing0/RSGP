function  [varargout] =covSEisoRew(varargin)
varargout = cell(max(1,nargout),1);
[varargout{:}] = covScaleRew({'covSE','iso',[]},varargin{:});

% global Reward_
% 
% if isnumeric(varargout{1}) 
%     K = varargout{1};
%     if size(K,1)==size(K,2) && size(K,1)==length(Reward_)      
%         indx = find(Reward_ ==0);  
%         MK = ones(size(K));
%         for k = 1:length(indx)
%                   kd = indx(k)+1; %min(indx(k)+1,size(K,1));
%                   if kd<size(K,2)
%                 MK(indx(k),kd:end) = 0;  MK(kd:end,indx(k)) = 0; 
%                   end
% 
%         end
%        MK = (MK+MK')./2 ;
%      varargout{1} = K.*MK;     %ishermitian(varargout{1})    
%    
%     % varargout{2} = varargout{2};
%     end
% 
% %      figure; imagesc(varargout{1});
% %      pause; 
% end
