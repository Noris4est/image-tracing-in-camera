%% Реализация двумерной свертки матрицы c переменным ядром
function [C] = special_conv2d(A,x_ref,y_ref,dxy,B,r_arr,add_value_n)
%dxy mkm - шаг дискретизации
%r_arr [mm, ... mm]
%xref yref mm - опорные координаты входного сегмента изображения A
%A - входная const матрица. C-возвращаемая матрица.
%B - матрица (переменное ядро)
%add_n - относительный параметр дополнения матрицы 0...1.
addValue=add_value_n*max(A,[],'all');%
%дополнение матрицы A для сохранения размерности при свертке. 
%Добавляются Nadd строк/столбцов со всех сторон.
Nadd_n =(size(B,1)-1)/2;
Nadd_m =(size(B,2)-1)/2;
[rows_A, columns_A]=size(A);
A_expansion=addValue*ones(rows_A+2*Nadd_n,columns_A+2*Nadd_m);%Матрица-заготовка
%Размещение матрицы А в центре матрицы-заготовки
A_expansion(Nadd_n+1:rows_A+Nadd_n,Nadd_m+1:columns_A+Nadd_m)=A;
Result=zeros(size(A));
dxy=dxy/1000;
Tr_y=dxy*(size(A,1)-1);
Tr_x=dxy*(size(A,2)-1);

x_cur=-Tr_x/2+x_ref;
y_cur=Tr_y/2+y_ref;
grid_n=1:size(B,1);
grid_m=1:size(B,2);
for i=1:rows_A
    for j=1:columns_A
        r_cur=sqrt(x_cur^2+y_cur^2);
        %{
        Интерполяция функции рассеяния точки ФРТ [point spread function PSF] 
        Исходный трехмерный массив содержит size(B,3) матрицы ФРТ.
        Каждая матрица соответствует виду ФРТ на определенном удалении от
        оптической оси (ОО) в плоскости изображения (квазифокальной плоскости
        объектива). Массив расстояний от ОО - r_arr
        Настоящая операция позволяет получить матрицу ФРТ для произвольной точки в
        плоскости изображения.
        %}
        B_cur=interp3(grid_m,grid_n,r_arr,B,grid_m,grid_n,r_cur);
        angle_cur=angle(x_cur+1i*y_cur);%определение углового полож. точки.
        if angle_cur<0
            angle_cur=angle_cur+2*pi;
        end
        angle_cur_deg=angle_cur*180/pi;
        %поворот матрицы ФРТ
        B_cur_rot=imrotate(B_cur,angle_cur_deg,'bilinear','crop');
        B_reverse1=flip(B_cur_rot,1);
        B_reverse=flip(B_reverse1,2);%Разворот матрицы по двум измерениям
        PartAe=A_expansion(i:i+size(B,1)-1,j:j+size(B,2)-1);
        Result(i,j)=sum(PartAe.*B_reverse,'all');
        x_cur=x_cur+dxy;
    end
    x_cur=-Tr_x/2+x_ref;
    y_cur=y_cur-dxy;
end
C=Result;
end