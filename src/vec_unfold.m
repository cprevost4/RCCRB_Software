function [J_12 J_13 J_23]   =     vec_unfold(dim)
%-------------------------------------------------------------------------%
% [J_12 J_13 J_23]=     vec_unfold(dim)
%
% This function generates the permutation matrix for passing from one vec-
% torization of a given unfolding to a vectorization of another unfolding.
%
% The inputs are:
% 
% - dim         : tensor dimensions.
%
% The outputs are:
% 
% - [J_12 J_13 J_23]        : permutation matrices.
%
% List of updates                 -     16/03/2015  -     R. C. Farias
%                                       Creation of the file 
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
ind =     zeros(dim);
l   =     0;
% Indexes in tensor format
for i     =     1:dim(1)
    for   j     =     1:dim(2)
        for     k     =     1:dim(3)
          l=l+2;
          ind(i,j,k)  =     l;
        end
    end
end

% Index for matricization
m_ind     =     prod(dim)./dim;

% First unfolding and its vectorization
ind_1     =    reshape(permute(ind,[1,2,3]),dim(1),m_ind(1));
v_1       =    ind_1(:);
% Second unfolding and its vectorization
ind_2     =    reshape(permute(ind,[2,1,3]),dim(2),m_ind(2));
v_2       =    ind_2(:);
% Third unfolding and its vectorization
ind_3     =    reshape(permute(ind,[3,1,2]),dim(3),m_ind(3));
v_3       =    ind_3(:);

% Mappings between vectorizations
[~,mapping_12] =    ismember(v_1,v_2);
[~,mapping_13] =    ismember(v_1,v_3);
[~,mapping_23] =    ismember(v_2,v_3);

% Initialization of the permutation matrices
J_12=     zeros(prod(dim),prod(dim));
J_13=     zeros(prod(dim),prod(dim));
J_23=     zeros(prod(dim),prod(dim));

% Filling with ones in the right positions
indexes   =     1:prod(dim);
i_12      =     sub2ind([prod(dim) prod(dim)],mapping_12,indexes');
i_13      =     sub2ind([prod(dim) prod(dim)],mapping_13,indexes');
i_23      =     sub2ind([prod(dim) prod(dim)],mapping_23,indexes');
J_12(i_12)=     1;
J_13(i_13)=     1;
J_23(i_23)=     1;
%-------------------------------------------------------------------------%
end