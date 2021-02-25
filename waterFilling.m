function lambdaF = waterFilling(CNR,deltaF,K0)
    if any(CNR == 0)
        lambda = 0;
    else
        lambda = ((1 + sum(deltaF./(K0*CNR)))*log(2))/(deltaF*length(CNR));
    end
    for i = 1:length(CNR)
        initP = (lambda*deltaF/log(2)) - deltaF/(K0*CNR(i));
        while any( initP < 0 )
            negIndex        = initP <= 0;
            posIndex        = initP >  0;
            NkRem           = nnz(posIndex);
            CNRRem          = CNR(posIndex);
            powAllcTemp = (lambda*deltaF/log(2)) - deltaF/(K0*CNR(i));
            initP(negIndex) = 0;
            initP(posIndex) = powAllcTemp;
        end
        P(i)              = initP;    
    end
    for i = 1:length(CNR)
        lambdaF(i) = ((P(i) + deltaF/(K0*CNR(i)))*log(2))/deltaF;
    end
    lambdaF = lambdaF(1,1);
end