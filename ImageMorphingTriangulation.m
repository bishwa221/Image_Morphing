function M = ImageMorphingTriangulation(warp_frac,dissolve_frac)

if nargin < 1
    warp_frac = .5;
end

if nargin < 2
    dissolve_frac= warp_frac; 
end


% ream images
I = im2double(imread('a.png'));
J = im2double(imread('c.png'));

% load mat file with points, variables Ip,Jp
 load('points.mat');
 
% initialize ouput image (morphed)
M = zeros(size(I));

%  Triangulation (on the mean shape)
% We create triangles here with an intermediate shape. MeanShape.

MeanShape = (1/2)*Ip+(1/2)*Jp;
TRI = delaunay(MeanShape(:,1),MeanShape(:,2));


% number of triangles
TriangleNum = size(TRI,1); 

% find coordinates in images I and J
CordInI = zeros(3,3,TriangleNum);
CordInJ = zeros(3,3,TriangleNum);
%Create duplicate triangles in TRI in Array CordInI and CordInJ
%Now image I and J have same set of delaunay triangles as TRI (intermediate image).
for i =1:TriangleNum
  for j=1:3
    
    CordInI(:,j,i) = [ Ip(TRI(i,j),:)'; 1];
    CordInJ(:,j,i) = [ Jp(TRI(i,j),:)'; 1]; 
    
  end
end

% create new intermediate shape according to warp_frac
Mp = (1-warp_frac)*Ip+warp_frac*Jp; 

 
% create a grid for the morphed image
[x,y] = meshgrid(1:size(M,2),1:size(M,1));

% for each element of the grid of the morphed image, find  which triangle it falls in
%TM = triangle ko index
TM = tsearchn([Mp(:,1) Mp(:,2)],TRI,[x(:) y(:)]);
TM(isnan(TM))=[];

t = 0;
% YOUR CODE STARTS HERE
for i = round(TM)
t = t + 1;    
vert_ind = TRI(i);%%%%make_it_t
cd = zeros(1,6);
for j = 1:3
cd(2*j-1) = Mp(vert_ind(j),1);
cd(2*j) = Mp(vert_ind(j),2);
end
x_123 = [cd(1),cd(3),cd(5)];
y_123 = [cd(2),cd(4),cd(6)];
A = [x_123;y_123;ones(1,3)];
pt = [x(t);y(t);1];
bay_cord = A\pt;
Ai = CordInI(:,:,i);
Aj = CordInJ(:,:,i);
homo_I = Ai*bay_cord;
cors_I = round(homo_I(1:2)/homo_I(3));
homo_J = Aj*bay_cord;
cors_J = round(homo_J(1:2)/homo_J(3));
         



IndI = sub2ind(size(I),cors_I(2),cors_I(1));
IndJ =  sub2ind(size(J),cors_J(2),cors_J(1));
IndM =  sub2ind(size(M),y(t),x(t));

% YOUR CODE ENDS HERE



% cross-dissolve
M(IndM)=uint8((1-dissolve_frac)* I(IndI)+ dissolve_frac * J(IndJ));
end


%figure(100);
%subplot(1,3,1);
%imshow(I);
%hold on;
%triplot(TRI,Ip(:,1),Ip(:,2))
%hold off;
%title('First')

%subplot(1,3,2);
%imshow(M);
%hold on;
%triplot(TRI,Jp(:,1),Jp(:,2))
%hold off
%title('Morphed')

%subplot(1,3,3);
%imshow(J);
%hold on;
%triplot(TRI,Jp(:,1),Jp(:,2))
%hold off
%title('Second')

end