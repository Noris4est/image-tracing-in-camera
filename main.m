clc
clear
close all
%% моделирование ФРТ объектива с характерной комой
%задаем ФРТ для двух точек в плоскости изоражения
%характеристики ФРТ - СКО продольная и поперечная радиальному направлению
sigma_across_arr=[1,1];%mkm
sigma_along_arr=[1,4];%mkm
r_arr=[0,0.07];%mm координаты точек
dxy=0.5;%мкм шаг сетки разбиения поверхности ФРТ
Tr_x=20;%мкм диапазон рассмотрения функции ФРТ по осям x и y
Tr_y=20;%мкм

gridx=-Tr_x/2:dxy:Tr_x/2;
gridy=-Tr_y/2:dxy:Tr_y/2;

[x,y]=meshgrid(gridx,gridy);%формирование координатной сетки
PSF_M=zeros(length(gridy),length(gridx),length(r_arr));
%формирование слоев матриц ФРТ
for i=1:length(r_arr)
    sigma_across=sigma_across_arr(i);
    sigma_along=sigma_along_arr(i);
    PSF_M(:,:,i)=exp(-(x.^2/sigma_along^2+y.^2/sigma_across^2));
    PSF_M(:,:,i)=PSF_M(:,:,i)/sum(PSF_M(:,:,i),'all');%нормировка sum=1
end
figure;
subplot(1,2,1);contourf(PSF_M(:,:,1));title('PSF contour in axis point');
subplot(1,2,2);contourf(PSF_M(:,:,2));title('PSF contour in field edge');

figure;mesh(x,y,PSF_M(:,:,1));title('PSF mesh in axis point');
figure;mesh(x,y,PSF_M(:,:,2));title('PSF mesh in field edge');
%% Формирование исходного изображения
N=128;%размеры в пикселах
M=128;
dn=25;%шаг в пикселах между центрами элементов изображения по измерениям
dm=25;
Im0=comb2d(N,M,dn,dm);%формирование матрицы - двумерной функции comb
figure;imagesc(Im0);title('2d comb function');
%размещение в ненулевых узлах матрицы rect-образных импульсов
Im0=conv2(Im0,ones(10),'same'); 

x_ref=0;%опорные координаты изображения Im0
y_ref=0;
%% свертка исх. изобр. с переменным ядром  
Im_out=special_conv2d(Im0,x_ref,y_ref,dxy,PSF_M,r_arr,0);
Tr_x_im=(M-1)*dxy/1000;
Tr_y_im=(N-1)*dxy/1000;
grid_im_x=-Tr_x_im/2+x_ref:dxy/1000:Tr_x_im/2+x_ref;
grid_im_y=-Tr_y_im/2+y_ref:dxy/1000:Tr_y_im/2+y_ref;

[x_im,y_im]=meshgrid(grid_im_x,grid_im_y);
figure;mesh(x_im,y_im,Im0);title('Image in f-plane ideal OS');
figure;mesh(x_im,y_im,Im_out);title('Image in f-plane real OS');
figure;
subplot(1,2,1);imagesc(Im0);title('Image in f-plane real OS');
subplot(1,2,2);imagesc(Im_out);title('Image in f-plane ideal OS');


