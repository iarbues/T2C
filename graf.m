function graf(DGP,P,R,b,tau,dist)

nombre=['outDGP' num2str(DGP) '/res_DGP' num2str(DGP) '_P' num2str(P) '_R' num2str(R) '_b' num2str(10*b) '_tau' num2str(tau) '_d' num2str(dist)]
load(nombre,'MSEmM_pval','MSEMm_pval','MSEnaive_0_pval','MSEnaive_1a_pval','MSEnaive_1b_pval');

[f1,pval1] = ecdf(MSEmM_pval);
selec=(pval1<=0.2);
plot(pval1(selec),f1(selec),':','LineWidth',3)
hold on;

[f2,pval2] = ecdf(MSEMm_pval); 
selec=(pval2<=0.2);
plot(pval2(selec),f2(selec),'r-.','LineWidth',3)

[f3,pval3] = ecdf(MSEnaive_0_pval);
selec=(pval3<=0.2);
plot(pval3(selec),f3(selec),'k','LineWidth',1)

[f4,pval4] = ecdf(MSEnaive_1a_pval);
selec=(pval4<=0.2);
plot(pval4(selec),f4(selec),'g--','LineWidth',1)

[f5,pval5] = ecdf(MSEnaive_1b_pval);
selec=(pval5<=0.2);
plot(pval5(selec),f5(selec),'y','LineWidth',3)

line([0 0.2],[0 0.2],'Color','k')
grid on;

legend('MSEmM','MSEMm','MSEnaive 0','MSEnaive 1a','MSEnaive 1b')

%title(['P=' num2str(P) ' R=' num2str(R) ' b=' num2str(b) ' tau=' num2str(tau) ' dist=' num2str(dist)]);

save([nombre '_datos']);
print('-deps',[nombre '.eps']);
print([nombre '.fig'])
print('-djpeg',[nombre '.jpg'])