
function graf(DGP,P,R,b,tau,dist,selmet,format,w,filename)

nombre=['outDGP' num2str(DGP) '/res_DGP' num2str(DGP) '_P' num2str(P) '_R' num2str(R) '_b' num2str(10*b) '_tau' num2str(tau) '_d' num2str(dist)]
load(nombre,'MSEmM_pval','MSEMm_pval','MSEnaive_0_pval','MSEnaive_1a_pval','MSEnaive_1b_pval','MSEnaive_2_pval','MSEnaive_3_pval');

if selmet(1)
    
    [f1,pval1] = ecdf(MSEmM_pval);
    selec=(pval1<=0.2);
    plot(pval1(selec),f1(selec),format{1},'LineWidth',w(1))
    hold on;
    
end

if selmet(2)
    
    [f2,pval2] = ecdf(MSEMm_pval);
    selec=(pval2<=0.2);
    plot(pval2(selec),f2(selec),format{2},'LineWidth',w(2))
    
end

if selmet(3)
    
    [f3,pval3] = ecdf(MSEnaive_0_pval);
    selec=(pval3<=0.2);
    plot(pval3(selec),f3(selec),format{3},'LineWidth',w(3))
    
end

if selmet(4)
    
    [f4,pval4] = ecdf(MSEnaive_1a_pval);
    selec=(pval4<=0.2);
    plot(pval4(selec),f4(selec),format{4},'LineWidth',w(4))
    
end

if selmet(5)
    
    [f5,pval5] = ecdf(MSEnaive_1b_pval);
    selec=(pval5<=0.2);
    plot(pval5(selec),f5(selec),format{5},'LineWidth',w(5))
    
end

if selmet(6)
    
    [f6,pval6] = ecdf(MSEnaive_2_pval);
    selec=(pval6<=0.2);
    plot(pval6(selec),f6(selec),format{6},'LineWidth',w(6))
    
end

if selmet(7)
    
    [f7,pval7] = ecdf(MSEnaive_3_pval);
    selec=(pval7<=0.2);
    plot(pval7(selec),f7(selec),format{7},'LineWidth',w(7))
    
end


line([0 0.2],[0 0.2],'Color','k')
grid on;

noms={'MSEmM','MSEMm','MSEnaive 0','MSEnaive 1a','MSEnaive 1b','MSEnaive 2','MSEnaive 3'};

legend(noms(selmet));

%title(['P=' num2str(P) ' R=' num2str(R) ' b=' num2str(b) ' tau=' num2str(tau) ' dist=' num2str(dist)]);

%save([nombre '_datos']);
print('-deps',[filename '.eps']);
print([filename '.fig'])
print('-djpeg',[filename '.jpg'])