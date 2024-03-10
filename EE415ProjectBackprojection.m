%Muhammed Saadeddin Koçak 2232346

clearvars -except projectionmatrix sizeofimage

[numberofbeams stepsize]=size(projectionmatrix);
imagematrix=zeros(sizeofimage);

%Loading Image
%Realize original image is loaded only for comparison. In backprojection
%code A matrix is never used.
myimage=load("lena.mat");
C=struct2cell(myimage);%Cell
A=cell2mat(C);%Matrix
figure
imagesc(A);
title ("Original Image");
colormap gray;
%Filtering
%Triangular Filter
% triangfilter=triang(numberofbeams);%Designing filter in Fourier Domain
% for o=1:stepsize
% filteredprojectionmatrix(:,o)=ifft2((triangfilter).*(fft2(projectionmatrix(:,o))));
% end
% projectionmatrix=filteredprojectionmatrix;
%Gaussian Filter
inpargu=1:1:numberofbeams;
gausfilter=normpdf(inpargu,numberofbeams/2,20);
gausfilter=gausfilter';
for o=1:stepsize
filteredprojectionmatrix(:,o)=ifft2((gausfilter).*(fft2(projectionmatrix(:,o))));
end
projectionmatrix=filteredprojectionmatrix;
%Plotting Filter
% figure
% plot(gausfilter);
%Defining t, teta, x, and y arrays according to beam number and step size
t=linspace(-sizeofimage*sqrt(2)/2,sizeofimage*sqrt(2)/2,numberofbeams);
teta=linspace(0,180-(180/stepsize),stepsize);
x=linspace(-sizeofimage/2,sizeofimage/2,sizeofimage+1);
y=linspace(-sizeofimage/2,sizeofimage/2,sizeofimage+1);
 for i=1:(length(teta))
     for j=1:(length(t))
%Clearing dynamic memory
relevantmatrix=[];
rowdata=[]; columndata=[]; midpointx=[];
midpointy=[]; distance=[];
%Finding intersection points
    yfromx=(t(j)/sind(teta(i)))-(x.*cotd(teta(i)));
    xfromy=(t(j)/cosd(teta(i)))-(x.*tand(teta(i)));
%Detecting relevant x and y coordinates
m=1;
for g=1:(length(x))
    if (((-sizeofimage/2)<=yfromx(g))&&(yfromx(g)<=sizeofimage/2))
        relevantmatrix(m,1)=x(1,g);
        relevantmatrix(m,2)=yfromx(1,g);
        m=m+1;
    else
    end
end
for k=1:(length(y))
    if (((-sizeofimage/2)<=xfromy(k))&&(xfromy(k)<=sizeofimage/2))
        relevantmatrix(m,1)=xfromy(1,k);
        relevantmatrix(m,2)=y(1,k);
        m=m+1;
    else
    end
end
uniquerelevantmatrix=unique(relevantmatrix,'rows');%Deleting same x and y points
sortedrelevantmatrix=sortrows(uniquerelevantmatrix);%Sorting relevant x and y points
sizerelevant=size(uniquerelevantmatrix,1);%Detecting if there are relevant points
if sizerelevant==0 %If no relevant points do nothing
else
[rowsizeofURM,columnsizeofURM]=size(uniquerelevantmatrix);
if (rowsizeofURM>1)%If there are more than one relevant point then calculate attenuation
%Calculating distance and midpoints to use them for detecting row and
%column data
for r=1:(rowsizeofURM-1)
    distance(r)=sqrt(((uniquerelevantmatrix(r,1)-uniquerelevantmatrix(r+1,1))^2)+((uniquerelevantmatrix(r,2)-uniquerelevantmatrix(r+1,2))^2));
    midpointx(r)=(uniquerelevantmatrix(r,1)+uniquerelevantmatrix(r+1,1))/2;
    midpointy(r)=(uniquerelevantmatrix(r,2)+uniquerelevantmatrix(r+1,2))/2;
end
%Finding corresponding column and row data
rowdata=(sizeofimage/2)-floor(midpointy);
columndata=(sizeofimage/2)+ceil(midpointx);
%Backprojecting one ray beam
sizeofdistancematrix=size(distance);
sizeofdistancearray=sizeofdistancematrix(1,2);
for h=1:sizeofdistancearray
    imagematrix(rowdata(1,h),columndata(1,h))=imagematrix(rowdata(1,h),columndata(1,h))+(projectionmatrix(j,i)*distance(1,h));
end
else%If there is only one relevant point beam is tangential, then do nothing
end
end
    end
 end
%Normalizing image matrix
imagematrix=real(imagematrix);
normimagematrix=(imagematrix-min(min(imagematrix)))/(max(max(imagematrix)));
%Displaying Image
figure
imagesc(normimagematrix);
title("Backprojected Image");
colormap gray;
%Quantitative Comparison
for i=1:sizeofimage
    for j=1:sizeofimage
        difference(i,j)=abs(normimagematrix(i,j)-A(i,j));
    end
end
%Displaying Differences Between Images
figure
imagesc(difference);
title("Difference Between Images");
colormap gray;
%Calculating Error
error=sum(sum(difference));
errorratio=error/(sizeofimage^2);