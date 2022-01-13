function [Im_out] = comb2d(rows,cols,Ti,Tj)
%функция формирует 2d распределение дельта-функций с заданным периодом
n0=fix((rows-1)/2)+1;
m0=fix((cols-1)/2)+1;
M_sample=zeros(rows,cols);
dm=0;%начальные значения абсолютного сдвига по направлениям
dn=0;
while (n0+dn<rows)&&(n0-dn>0)
    while (m0+dm<cols)&&(m0-dm>0)
        M_sample(n0+dn,m0+dm)=1;
        M_sample(n0-dn,m0+dm)=1;
        M_sample(n0+dn,m0-dm)=1;
        M_sample(n0-dn,m0-dm)=1;
        dm=dm+Tj;
    end
    dm=0;
    dn=dn+Ti;
end
    Im_out=M_sample;
end

