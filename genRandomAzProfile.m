function [azProfile, H_TX, H_RX, physH] = genRandomAzProfile(prm, ...
                                                            thetaMin, thetaMax, ...
                                                            BsSteer, RxSteer)

   % Generate azimuth profile demonstrating incoherent DFT dictionaries
   % BsArrayPos, RxArrayPos are NumBsElements x D, NumRxElements x D
   % D the dimension of the array
   
   azBins = prm.K;
   nScat = 2;

   thetaBins = thetaMin:(thetaMax-thetaMin)/(azBins-1):thetaMax;

   azProfile = zeros(azBins, 1);
   scatAngs = randperm(azBins, nScat);
    
   azProfile(scatAngs) = complex(1, 1) / sqrt(2);
   %H_TX is NumBsElements x K
   %H_RX is NumRxElements x K

%    H_TX = generateDFT(prm.NumBsElements, azBins + 1);
%    H_RX = generateDFT(prm.NumRxElements, azBins + 1);

   H_TX = zeros(prm.NumBsElements, azBins);
   H_RX = zeros(prm.NumRxElements, azBins);

   for k = 1:azBins
        H_TX(:, k) = BsSteer(prm.CenterFreq, thetaBins(k));
        H_RX(:, k) = RxSteer(prm.CenterFreq, thetaBins(k));
   end
   physH = H_RX * diag(azProfile) * H_TX.';

end
