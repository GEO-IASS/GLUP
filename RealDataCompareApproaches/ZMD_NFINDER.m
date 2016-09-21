function [E,E_indicie] = ZMD_NFINDER(R,p)
%[E,E_indicie] = ZMD_NFINDER(R,p)
%N_FINDER ALGORITHME(Endmember Extraction) 
%The main idea of N-FINDR is to assume that a (p - 1) dimensional volume formed by a simplex 
%with p vertices that are specified by the purest pixels is always larger than that formed by 
%another combination of p pixels.
%------------------------------------------------------
%Usage
%Inputs
%  R = 2D matrix of HSI data (L x N )
%  p = Number of endmembers 
%ourputs
%  E = Matrix of endmember (L x p )
%  E_indicie = Indicies of pure pixels in E
%
%  References:
%     [1]M. E. Winter, “N-finder: An algorithm for fast autonomous spectral endmember determination 
%        in hyperspectral data,” in Proc. Image Spectrom. V, 1999, vol. 3753, pp. 266–277.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parametres
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if (nargin == 0)
           error('Insufficient parameters');
       end
       if (nargin < 2)
           p = ZMD_VDhfc(R,10^(-5));
       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initializations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if isempty (R)
          error ('there is no data');
       else
          [L,N]= size (R);%L number of bands (channels)
                    %N number of pixels
       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Dimensionality reduction by PCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       r_bar = mean(R,2);
       R_zeromean = R - repmat (r_bar,1,N);
       [COEFF, ZSCORES] = princomp(R_zeromean.');
       R_pca = (ZSCORES(:,1:p-1))';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       E = zeros (L,p);
       E_indicie = zeros(1,N);
       E_init = zeros(p);
       E_init(1,:) = 1;
       IDX = randperm(N);
       IDX = IDX(1:p);
       E_init(2:p,:) = R_pca(:,IDX);
       Vol_init = abs(det(E_init))/(factorial(p-1));%Volume calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Stopping rule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       for k=1:p
           for i=1:N
               e_k = E_init (2:p,k);  
               E_init(2:p,k) = R_pca(:,i);
               Vol_new =  abs(det(E_init))/(factorial(p-1));
               if Vol_new > Vol_init
                  Vol_init = Vol_new;
                  IDX(k) = i;
               else 
                  E_init(2:p,k) = e_k;
               end
           end
       end

       E = R(:,IDX);
       E_indicie = IDX; 

       return







