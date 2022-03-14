function vf_mat = viewfactorBo( c ,outlayer)
%% viewfactorBo calculate the area * view factor eachother
%   c is the coordinate of the panel start and end points
%   [x1st, y1st, x1end, y1end; 
%    x2st, y2st, x2end, y2end;
%    ...]  with Clockwise direction
%% 
n = length(c);       % number of points (or lines)
%m = n*(n-1)/2;                 % total number of viewfactors

vf1 = zeros(n);
 h = 1;
for i = 1:(n-1)
       for j = 1:(n-i)
        A = c(i,(1:2));
        B = c(i,(3:4));
        C = c(i+j,(1:2));
        D = c(i+j,(3:4));
        AC = sqrt(sum((A-C).^2));
        AD = sqrt(sum((A-D).^2));
        BC = sqrt(sum((B-C).^2));
        BD = sqrt(sum((B-D).^2));
        if AC==0 && BD==0
            vf1(j+h,i) = 1;
        else
            vf1(j+h,i) = (AC+BD-BC-AD)/2;   % crossed string method  
        end
       end
     h = h+1;
end
%vf1(vf1<0)=0;         % transfer the negative number to zero
vf_mat1 = vf1+vf1';    % matrix of view factor
vf_mat1=vf_mat1-diag(diag(vf_mat1));        % the diag is 0

b = zeros(n);                          % change the external layer to 1 for sky
b(:,n) = outlayer;
b(n,:) = outlayer';

vf_mat = (1-b).*vf_mat1 + b;

%vf = nonzeros(vf1);   % list the none zero number of vf1
    
end

