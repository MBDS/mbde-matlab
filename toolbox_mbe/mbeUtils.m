classdef mbeUtils
    % A "namespace" for miscelaneous auxiliary static methods
    
	% -----------------------------------------------------------------------------
	% This file is part of MBDE-MATLAB.  See: https://github.com/MBDS/mbde-matlab
	% 
	%     MBDE-MATLAB is free software: you can redistribute it and/or modify
	%     it under the terms of the GNU General Public License as published by
	%     the Free Software Foundation, either version 3 of the License, or
	%     (at your option) any later version.
	% 
	%     MBDE-MATLAB is distributed in the hope that it will be useful,
	%     but WITHOUT ANY WARRANTY; without even the implied warranty of
	%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	%     GNU General Public License for more details.
	% 
	%     You should have received a copy of the GNU General Public License
	%     along with MBDE-MATLAB.  If not, see <http://www.gnu.org/licenses/>.
	% -----------------------------------------------------------------------------

    methods(Static,Access=public)
        %---------------------------------------------------------------
        function x=qchisq(P,n)
            % QCHISQ(P,N) - quantile of the chi-square distribution.
            if nargin<2
              n=1;
            end

            s0 = P==0;
            s1 = P==1;
            s = P>0 & P<1;
            x = 0.5*ones(size(P));
            x(s0) = -inf;
            x(s1) = inf;
            x(~(s0|s1|s))=nan;

            for ii=1:14
              dx = -(mbeUtils.pchisq(x(s),n)-P(s))./mbeUtils.dchisq(x(s),n);
              x(s) = x(s)+dx;
              if all(abs(dx) < 1e-6)
                break;
              end
            end
        end

        %---------------------------------------------------------------
        function F=pchisq(x,n)
            % PCHISQ(X,N) - Probability function of the chi-square distribution.
            if nargin<2
              n=1;
            end
            F=zeros(size(x));

            if rem(n,2) == 0
              s = x>0;
              k = 0;
              for jj = 0:n/2-1;
                k = k + (x(s)/2).^jj/factorial(jj);
              end
              F(s) = 1-exp(-x(s)/2).*k;
            else
              for ii=1:numel(x)
                if x(ii) > 0
                  F(ii) = quadl(@mbeUtils.dchisq,0,x(ii),1e-6,0,n);
                else
                  F(ii) = 0;
                end
              end
            end
        end % pchisq()
        
        function f=dchisq(x,n)
            % DCHISQ(X,N) - Density function of the chi-square distribution.
            if nargin<2
              n=1;
            end
            f=zeros(size(x));
            s = x>=0;
            f(s) = x(s).^(n/2-1).*exp(-x(s)/2)./(2^(n/2)*gamma(n/2));
        end % dchisq ()
        
        % -------------------
        function [S] = sigma_sampling_generate(MEAN,COV)
            % Draw deterministic samples (for UKF methods and alike)
            %
            
            if (any(eig(COV)<=1e-12))
                error('Ey! Non pos. def. covariance matrix.');
            end
            P_chol = chol(COV,'lower');
            L = length(MEAN);
            S = zeros(1+2*L, L);

            % make sure MEAN is a row vector:
            if (size(MEAN,1)>1)
                MEAN=MEAN';
            end

            % Typical values:
            alpha = 1e-3;
            k = 0;

            lambda = alpha^2*(L+k)-L;
            gamma = sqrt(L+ lambda);

            % 0'th:
            S(1,:) = MEAN;

            % 1:2L+1
            for i=1:L,
                S(2*i+0,:) = MEAN + gamma * P_chol(:,i)';
                S(2*i+1,:) = MEAN - gamma * P_chol(:,i)';
            end
        end % sigma_sampling_generate
        
        
        function [X,P, Pzx] = sigma_sampling_mean_and_cov(Xs, Zs,Zm)
            % Recover the mean and covariance matrix from a set of samples 
            % (each row in "Xs" is a sample). Optionally, if Z is present (Zs samples,
            % Zm mean), this also computes the cross covariance Cov(Z,X)
            Nsp = size(Xs,1); % 2L+1
            L = (Nsp-1)/2;

            nX = size(Xs,2);

            % Typical params:
            alpha = 1e-3;
            k = 0;
            beta = 2; % Optimal for Gaussians
            lambda = alpha^2*(L+k)-L;


            % Weights:
            wm = zeros(Nsp,1);  % mean weights
            wc = zeros(Nsp,1);  % cov weights

            wm(1) = lambda / (L+lambda);
            wc(1) = wm(1) + (1-alpha^2+beta);

            wm(2:Nsp) = ones(Nsp-1,1) * 0.5/(L+lambda);
            wc(2:Nsp) = wm(2:Nsp);

            % Mean:
            X = wm' * Xs;

            % Cov:
            P=zeros(nX,nX);
            for i=1:Nsp,
                Ax =  (Xs(i,:)-X);
                P = P + wc(i)* Ax'*Ax;
            end

            % Cross cov(z,x)
            if nargin > 1
                assert(nargin==3);
                nZ = length(Zm);
                Pzx=zeros(nZ,nX);
                for i=1:Nsp,
                    Az =  (Zs(i,:)-Zm);
                    Ax =  (Xs(i,:)-X);
                    Pzx = Pzx + wc(i)* Az'*Ax;
                end
            end

        end % sigma_sampling_mean_and_cov()
        
        
    end % end methods()
    
end

