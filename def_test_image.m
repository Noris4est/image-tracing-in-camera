function [TestIm] = def_test_image(N,M,dn,dm,base_n,base_m,type)
Im_sample=comb2d(N,M,dn,dm);%двумерный периодич. след. дельта функц.
switch type
    case 'rect'
        core=ones(base_n,base_m);
    case 'pyr2d'
        [m,n]=meshgrid(1:base_m,1:base_n);
        profile_m=lambda(2*(m-base_m/2)/base_m);
        profile_n=lambda(2*(n-base_n/2)/base_n);
        core=min(profile_m,profile_n);
    case 'pyr1d'
        [m,n]=meshgrid(1:base_m,1:base_n);
        core=lambda(2*(n-base_n/2)/base_n);
    case 'circ'
        base=max(base_n,base_m);
        [m,n]=meshgrid(1:base,1:base);
        nm0=base/2;
        r=sqrt((n-nm0).^2+(m-nm0).^2);
        core=circ(2*r/base);
    case 'sphere'
        base=max(base_n,base_m);
        [m,n]=meshgrid(1:base,1:base);
        nm0=base/2;
        r=sqrt((n-nm0).^2+(m-nm0).^2);
        r0=base/2;
        core=double(r<=r0).*sqrt(r0^2-r.^2);
end
TestIm=conv2(Im_sample,core,'same');
end

