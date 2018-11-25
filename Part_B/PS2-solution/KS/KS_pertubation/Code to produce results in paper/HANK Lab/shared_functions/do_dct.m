function [index_reduced,DC1,DC2,DC3] = do_dct(obj,mpar, level)
% This function performs a discrete cosine transformation of a provided
% vector or higher-dimensional object

% It calculates the number of coefficients that are necessary to describe
% the input and collapses it to that number.

% The index, containing the position of information and zeros in the
% original vector is stored. The later can be used to find the original
% input using the inverse of the discrete cosine transformation.


obj = reshape(obj,[mpar.nm, mpar.nk, mpar.nh]);
[X1,DC1] = mydct(obj,1); % do dct-transformation
[X2,DC2] = mydct(X1,2); % do dct-transformation
[X3,DC3] = mydct(X2,3); % do dct-transformation
% Determine the number of data points sufficient to describe the control
% space
XX=X3(:);
[~,ind] = sort(abs(XX),'descend');
i = 1;
while norm(XX(ind(1:i)))/norm(XX) < level
    i = i + 1;
end
needed = i;
% fprintf('%3.1f%% of the coefficients are sufficient to describe the control space\n',needed/numel(X1)*100)

% Set irrelevant information to 0
index_reduced = sort(ind(1:i)); 

