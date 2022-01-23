clc
clear
close all
%% Определение ФРТ объектива. Графическое отображения слоев ФРТ.
dxy=0.5;%мкм - шаг между отсчетами в матрице ФРТ
%Функция формирования пользовательской синтетический ФРТ 
[PSF_M,r_arr,x_psf,y_psf]=define_user_PSF(dxy);
figure;
subplot(1,2,1);contourf(PSF_M(:,:,1));title('PSF contour in axis point');
subplot(1,2,2);contourf(PSF_M(:,:,2));title('PSF contour in field edge');

figure;mesh(x_psf,y_psf,PSF_M(:,:,1));title('PSF mesh in axis point');
figure;mesh(x_psf,y_psf,PSF_M(:,:,2));title('PSF mesh in field edge');
%% Формирование исходного изображения
N=128;%размеры в пикселах
M=128;
dn=40;%период следования дельта-импульсов на исх. изоражении по осям
dm=40;
base_n=25;%размеры основания периодически повторяющегося элемента
base_m=25;
figure;imagesc(comb2d(N,M,dn,dm));title('2d comb function');
Im0=def_test_image(N,M,dn,dm,base_n,base_m,'pyr2d');
%{ 
1. Полагаем матрицу Im0 исходным изображением при учете, 
что шаг между элементами dxy по осям.
2. Половина диагонали изображения не должна выходить за границы области
определения ФРТ по координате поля. Максимальная координата поля 
сформированной ФРТ - максимальное (последнее) значение массива r_arr
r_max>sqrt(M^2+N^2)*dxy/2 
%}

%{
Изображение объекта может иметь достаточно большую размерость. Операция
свертки с переменным ядром будет занимать длительное время. Для удобста и
быстроты анализа формирования объективом изображения в некоторой области
пл. изображения можно использовать в качестве исх. изобр. некий малый 
сегмент, для которого заданы опорные координаты.
(это координаты центра сегмента в плоскости изобр относительно опт. оси ОС). 
%}
x_ref=0;%опорные координаты изображения Im0
y_ref=0;
%% свертка исх. изобр. с переменным ядром  
Im_ROS=special_conv2d(Im0,x_ref,y_ref,dxy,PSF_M,r_arr,0);
Tr_x_im=(M-1)*dxy/1000;
Tr_y_im=(N-1)*dxy/1000;
grid_im_x=-Tr_x_im/2+x_ref:dxy/1000:Tr_x_im/2+x_ref;
grid_im_y=-Tr_y_im/2+y_ref:dxy/1000:Tr_y_im/2+y_ref;

[x_im,y_im]=meshgrid(grid_im_x,grid_im_y);
figure;mesh(x_im,y_im,Im0);title('Image in f-plane ideal OS');
figure;mesh(x_im,y_im,Im_ROS);title('Image in f-plane real OS');
figure;
subplot(1,2,1);imagesc(abs(fftshift(fft2(Im0))));title('Spectrum image in f-plane ideal OS');
subplot(1,2,2);imagesc(abs(fftshift(fft2(Im_ROS))));title('Spectrum image in f-plane real OS');

figure;
subplot(1,2,1);imagesc(Im0);title('Image in f-plane real OS');
subplot(1,2,2);imagesc(Im_ROS);title('Image in f-plane ideal OS');
%% регистрация изображения матричным приемником излучения
%{
Объектив камеры формирует распределение облученности E(x,y) в фокальной
плоскости. Распределение облученности представляет собой изображение,
которое нужно зарегистрировать. Для регистрации в настоящем случае
используется матричный приемник излучения (МПИ). МПИ последовательно
выполняет следующие процедуры. 
1) Усреднение потока
2) Дискретизация с квантованием
С точностью до коэффициента можно первую процедуру выполнить путем свертки
с функцией двумерного прямоугольного импульса (rect). Ширина импульса будет
задаваться равной ширине пиксела МПИ. 
%}
a_elem=2;%мкм квадратный пиксел
T_elem=a_elem;%fill factor 100%
M_av_im_sens=ones(fix(a_elem/dxy));%матрица усреднения
Im_ROS_av=conv2(Im_ROS,M_av_im_sens,'same');
%выборка усредненных потоков
U_MPI=Im_ROS_av.*comb2d(128,128,T_elem/dxy,T_elem/dxy);
figure;imagesc(U_MPI);title('image sensor signal');
%определение спектра сигнала с МПИ
F_U_MPI=fft2(U_MPI);
figure;imagesc(abs(fftshift(F_U_MPI)));title('Spectrum signal from im. sensor');






