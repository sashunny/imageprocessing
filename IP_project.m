%importing the image
image=imread('/home/ajsharma/Pictures/IP.jpg');
img=rgb2gray(image);
id=im2double(img);
subplot(3,4,1);imshow(img);title('Original Image');

%addition
add = img+30; 
subplot(3,4,2);imshow(add);title('Addition');

%subtraction
sub = img-15;
subplot(3,4,3);imshow(sub);title('Subtraction');

%multiplication
mul = img.*2;
subplot(3,4,4);imshow(mul);title('Multiplication');

%division
div = img./2;
subplot(3,4,5);imshow(div);title('Division');

%contrast stretching
%stretching the pixel values from 3 to 155 to the range of 0 to 255
con_stretch=imadjust(img,[0.01 0.6],[0 1]);
subplot(3,4,6);imshow(con_stretch);title('Cont. Stretch');

%thresholding
level=graythresh(img); %selects the best value for thresholding
i_thresh = im2bw(img,level); %using the im2bw function for thresholding
subplot(3,4,7);imshow(i_thresh);title('Thresholding');

%image negative
i_neg=255-img;
subplot(3,4,8);imshow(i_neg);title('Negative');

%log and power law transformation
i_log=id;
i_power=id;
[r c]=size(id);
c1=1.3; %log transfor constant
c2=1.5; %power law transform constant
power=2; 
for i = 1:r
    for j = 1:c
        i_log(i,j)=c1*log(1+id(i,j));
        i_power(i,j)=c2*id(i,j)^power;
    end
end
subplot(3,4,9);imshow(i_log);title('Log Transform');
subplot(3,4,10);imshow(i_power);title('Power Transform');

%gray level slicing
i_slicea=img;
i_sliceb=img;
[row1, col1]=size(img);
for i = 1:row1
    for j = 1:col1
       if i_slicea(i,j)>=150 && i_slicea(i,j)<230
           i_slicea(i,j)=255;
           i_sliceb(i,j)=255;
       else
           i_slicea(i,j)=0;
       end
    end
end
subplot(3,4,11);imshow(i_slicea, []);title('Gray Slicing A');
subplot(3,4,12);imshow(i_sliceb);title('Gray Slicing B');

%bit plane slicing
[pl1,pl2,pl3,pl4,pl5,pl6,pl7,pl8]=bitplane_code(img);
figure;
    subplot(3,3,1);imshow(pl1);title('1st pln');
    subplot(3,3,2);imshow(pl2);title('2nd pln');
    subplot(3,3,3);imshow(pl3);title('3rd pln');
    subplot(3,3,4);imshow(pl4);title('4th pln');
    subplot(3,3,5);imshow(pl5);title('5th pln');
    subplot(3,3,6);imshow(pl6);title('6th pln');
    subplot(3,3,7);imshow(pl7);title('7th pln');
    subplot(3,3,8);imshow(pl8);title('8th pln');
    rec=pl1+pl2*2+pl3*4+pl4*8+pl5*16+pl6*32+pl7*64+pl8*128;
    subplot(3,3,9);imshow(uint8(rec));title('Original Img')

%histogram equalization
hist=imhist(img);
i_eq=histeq(img);
hist_eq=imhist(i_eq);
figure;
    subplot(2,2,1);imshow(img);title('Original Image');
    subplot(2,2,2);imshow(hist);title('Original Histogram');
    imhist(img)
    subplot(2,2,3);imshow(i_eq);title('Equalized Image');
    subplot(2,2,4);imshow(hist_eq);title('Equalized Histogram');
    imhist(i_eq)

    
% STATISTICAL NON-LINEAR FILTERS
%minimum filter
%maximum filter
%median filter
B=zeros(size(img));
C=zeros(size(img));
D=zeros(size(img));
modifyA=padarray(img,[1 1]);
        x=[1:3]';
        y=[1:3]';
       
for i= 1:size(modifyA,1)-2
    for j=1:size(modifyA,2)-2 
       window=reshape(modifyA(i+x-1,j+y-1),[],1);
       B(i,j)=min(window);
       C(i,j)=max(window);
       D(i,j)=median(window);
    end
end

B=uint8(B);
C=uint8(C);
D=uint8(D);
figure;
    subplot(2,2,1);imshow(img),title('ORIGINAL IMAGE');
    subplot(2,2,2);imshow(B),title('IMAGE AFTER MIN FILTERING');
    subplot(2,2,3);imshow(C),title('IMAGE AFTER MAX FILTERING');
    subplot(2,2,4);imshow(D),title('IMAGE AFTER MEDIAN FILTERING');

%SMOOTHENING FILTER
%Mean Filter
mf=ones(3,3)/9;
i_smooth=filter2(mf,id);
figure;
    subplot(1,2,1);imshow(id),title('Original Image');
    subplot(1,2,2);imshow(i_smooth),title('After Smoothening');

%SHARPENING FILTER
%Laplacian Filter
h=fspecial('laplacian',0.2);
sharp=filter2(h,id);
i_sharp=id-sharp;
figure;
    subplot(1,2,1);imshow(id),title('Original Image');
    subplot(1,2,2);imshow(i_sharp),title('After Sharpening');
    
%Sobel Filter
h_sobel=fspecial('sobel');
i_sobel=filter2(h_sobel,id);
figure;
    subplot(1,2,1);imshow(id),title('Original Image');
    subplot(1,2,2);imshow(i_sobel),title('Soble Filter applied');
    
%Prewitt filter
h_prewitt=fspecial('prewitt');
i_prewitt=filter2(h_prewitt,id);
figure;
    subplot(1,2,1);imshow(id),title('Original Image');
    subplot(1,2,2);imshow(i_prewitt),title('Prewitt Filter applied');

%Robert Cross filter
i_robert=edge(img,'Roberts');
figure;
    subplot(1,2,1);imshow(img),title('Original Image');
    subplot(1,2,2);imshow(i_robert),title('Roberts Filter applied');
    
% ### MORPHOLOGICAL OPERATOR ### %
imgb=imbinarize(imread('/home/ajsharma/Pictures/IPB.jpg'));
%Erosion
se=[0 1 0;
    1 1 1;
    0 1 0];
i_erode=imerode(imgb,se);
%Dilation
i_dil=imdilate(imgb,se);
%edge extraction
i_ex=i_dil-imgb;
figure;
    subplot(2,2,1);imshow(imgb);title('Original Image');
    subplot(2,2,2);imshow(i_erode);title('Erosion');
    subplot(2,2,3);imshow(i_dil);title('Dilation');
    subplot(2,2,4);imshow(i_ex);title('Edge Extracted');
    
 %opening
 i_open=imdilate(imerode(imgb,se),se);
 %closing
 i_close=imerode(imdilate(imgb,se),se);
 figure;
    subplot(1,3,1);imshow(imgb);title('Original Image');
    subplot(1,3,2);imshow(i_open);title('Opening');
    subplot(1,3,3);imshow(i_close);title('Closing');
    
 %hit and miss transform
 i_hitmiss=bwhitmiss(imgb,se);
 
 %thinning
 i_thin=imgb-i_hitmiss;
 
 %thickening
 i_thick=imgb+i_hitmiss;
 figure;
    subplot(2,2,1);imshow(imgb);title('Original Image');
    subplot(2,2,2);imshow(i_hitmiss);title('Hit Miss Transform');
    subplot(2,2,3);imshow(i_thin);title('Thinning');
    subplot(2,2,4);imshow(i_thick);title('Thickening');

    
 %#### IMAGE TRANSFORMS ####%
 %DFT of an image
 fourier=fft2(img);
 %DCT of an image
 cosine=dct2(img);
 figure;
    subplot(1,3,1);imshow(img);title('Original Image');
    subplot(1,3,2);imshow(real(fourier));title('Fourier Transform');
    subplot(1,3,3);imshow(real(cosine));title('Cosine Transform');
    
    
    
%###FUCTION FOR BIT PLANE SLICING####%   
function [pl1 pl2 pl3 pl4 pl5 pl6 pl7 pl8]=bitplane_code(img)
img=rgb2gray(imread('/home/ajsharma/Pictures/IP.jpg'));
[row col]=size(img);
b=zeros(row,col,8);

for k=1:8
    for i=1:row
        for j=1:col
            b(i,j,k)=bitget(img(i,j),k);
        end
    end
end
pl1=b(:,:,1);
pl2=b(:,:,2);
pl3=b(:,:,3);
pl4=b(:,:,4);
pl5=b(:,:,5);
pl6=b(:,:,6);
pl7=b(:,:,7);
pl8=b(:,:,8);
end