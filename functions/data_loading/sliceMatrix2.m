function data=sliceMatrix2(M,slicewidth)
% M - nxm data matrix
% slicewidth - number of matrix columns to group together as a slice
% data - (n*slicewidth)xm matrix consisting of all the slices
% v.2 puts the delays AFTER the current time instead of before. 

[n,m]=size(M);

data=zeros(n*slicewidth,m);

% pad BEGINNING with zeros (so first slice is mostly zeros plus the first
% actual data point)
M=[zeros(n,slicewidth-1),M];

for i=1:m
    slice=M(:,i:i+slicewidth-1);
    data(:,i)=slice(:);
end

end