function bnd=detbndy(in)
%DETBNDY compute a boundary segment list for a FEM domain
% DETBNDY bnd=detbndy(e);
%         This function computes a boundary for the FEM domain
%         described a file containing element connectivity list (e).
%         It uses sparse matrix techniques to determine the element
%         edges on the boundary of the FEM domain.
%
% Input:  ele -  element list; 3 (.tri) or 4 (.ele) columns wide
% Output: bnd -  a 2-column list of boundary-node numbers, returned
%                to the local workspace
%
%         The output boundary list are pairs of node numbers, not 
%         coordinates, describing the edges of elements on the 
%         exterior of the domain, including islands.  The segments 
%         are not connected.
%
%         Call as: bnd=detbndy(e);
%
% Written by : Brian O. Blanton at The University of North Carolina 
%              at Chapel Hill, Mar 1995.
%

% DEFINE ERROR STRINGS
err1='Only one input argument to DETBNDY. Type "help detbndy"';
err2='Element list passed to DETBNDY does not have 3 or 4 columns';

% check argument list
if nargin~=1
   error(err1);
end

%if nargout~=1
%   disp(err3);
%   return
%end

% Check size of element list
[~,ncol]=size(in);
if ncol < 3 || ncol > 4
   error(err2);
elseif ncol==4
   in=in(:,2:4);
end

% Form (i,j) connection list from .ele element list
%
%i=[in(:,1);in(:,2);in(:,3)];
%j=[in(:,2);in(:,3);in(:,1)];
i=in;
j=circshift(in,[0 -1]);

% this sequence will get a sorted, ordered edge list, [a jjj] 
%[a,b]=sort([i j]');
%ii=a(1,:)';
%jj=a(2,:)';
%[a,b]=sort(ii);
%jjj=jj(b);
%[a jjj]

% Form the sparse adjacency matrix and add transpose.
%
n = max(max(i(:)),max(j(:)));
A = sparse(i,j,1,n,n);
A = A + A';

% Consider only the upper part of A, since A is symmetric
% 
%A=A.*triu(A);
A=triu(A);

% The boundary segments are A's with value == 1
% Interior segments (shared by 2 elements) are at value == 2
%B=A==1;
%Bi=A==2;

% Extract the row,col from B for the boundary list.
%
[ib,jb,s]=find(A==1);
bnd=[ib(:),jb(:)];

%[ibi,jbi,si]=find(A==2);
%intsegs=[ibi(:),jbi(:)];

% Output .bnd list
%
%res=input('Output .bnd boundary list (y/n) ? [n] ','s');
%if strcmp(res,'y')
%   fname=input('Enter .bnd filename [bnd.dat] ','s');
%   if isempty(fname),fname='bnd.dat';,end
%   fp=fopen(fname,'w');
%   for i=1:length(bnd)
%      fprintf(fp,'%d %d\n', bnd(i,1),bnd(i,2));  
%   end
%end
%
%LabSig  Brian O. Blanton
%        Department of Marine Sciences
%        12-7 Venable Hall
%        CB# 3300
%        University of North Carolina
%        Chapel Hill, NC
%                 27599-3300
%
%        brian_blanton@unc.edu
%

