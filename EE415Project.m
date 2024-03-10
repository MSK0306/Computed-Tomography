%Muhammed Saadeddin Koçak 2232346

clear
%Loading Image
myimage=load("square.mat");
C=struct2cell(myimage);%Cell
A=cell2mat(C);%Matrix
% A=[1 2;3 4];
%Getting Size
sizeofmatrix=size(A);
sizeofimage=sizeofmatrix(1,1);%Assuming image is square
%Getting beam number and step size
numberofbeams=300;%Number of beams
stepsize=180;%Number of steps
%Defining t, teta, x, and y arrays according to beam number and step size
t=linspace(-sizeofimage*sqrt(2)/2,sizeofimage*sqrt(2)/2,numberofbeams);
teta=linspace(0,180-(180/stepsize),stepsize);
x=linspace(-sizeofimage/2,sizeofimage/2,sizeofimage+1);
y=linspace(-sizeofimage/2,sizeofimage/2,sizeofimage+1);
%Taking projections
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
if sizerelevant==0 %If no relevant points exist return total attenuation
    totalattenuation=0;
    projectionmatrix(j,i)=totalattenuation;
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
%Summing all attenuation
sizeofdistancematrix=size(distance);
sizeofdistancearray=sizeofdistancematrix(1,2);
totalattenuation=0;%Clearing dynamic data from previous calculation
for h=1:sizeofdistancearray
    totalattenuation=totalattenuation+(A(rowdata(1,h),columndata(1,h))*distance(1,h));
end
projectionmatrix(j,i)=totalattenuation;%Assigning attenuation to matrix
else%If there is only one relevant point beam is tangential, then attenuation is zero
    totalattenuation=0;
    projectionmatrix(j,i)=totalattenuation;
end
end
    end
 end
 %Projection at a single angle
 figure
 plot(projectionmatrix(:,45));