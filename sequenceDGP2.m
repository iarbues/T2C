
dist=1; % distribution for the bootstrap



for P=80 % OOS length

    for tau=1               % horizon
        
        if tau==1,
            maxb=0.4;
        else
            maxb=0.8;
        end

        for b=[maxb 0]          % =0, null


            R=P;                % equal length
            T=P+R;
            lambda=P/T;

            modelsDGP2;         % models of both classes
            
            pval_comp;          % program that gets the p-values

            nombre0='outDGP2/res_DGP2_P'

            nombre=[nombre0 num2str(P) '_R' num2str(R) '_b' num2str(10*b) '_tau' num2str(tau) '_d' num2str(dist)];

            save(nombre);

        end
    end

end

