function [signalFiltered2D,Sigm]= funcDRT(signalSim2D,nRotate,Dim)
[M,N] = size(signalSim2D);
signalRot2D = zeros(2*M,2*N,nRotate);
signalRevRot2D = 0;
sigmaThr = 1;
alpha = 1;
for k = 1:nRotate
         temp =(imrotate(signalSim2D,360*(k-1)/nRotate));
         [Mtemp,Ntemp] = size(temp);
         padM = 2*M - Mtemp;
         padN = 2*N - Ntemp;
         signalRot2D = padarray(temp,[padM padN],'post');
         [U,Sigm,V] = svd(signalRot2D);
         signalEst2D = 0;
       %  rateSigm = abs(diff(diag(Sigm)));
       %  idx = find(rateSigm<sigmaThr);
         %P = ceil(1*idx(1));
         for p = 1:Dim
          signalEst2D =   signalEst2D  + (Sigm(p,p)^alpha)*U(:,p)*V(:,p)';
         end
         signalEst2D = (Sigm(1,1)^(1-alpha))*signalEst2D;
         diagSigm = diag(Sigm);
         diagSigmVec(k) = sum(diagSigm(1:Dim));
         signalRot2D = signalEst2D(1:end-padM, 1:end-padN);
         temp = imrotate(signalRot2D,-360*(k-1)/nRotate);
         [Mtemp,Ntemp] = size(temp);
         cropM = floor(0.5*(Mtemp - M));cropN = floor(0.5*(Ntemp - N));
         signalRevRot2D = signalRevRot2D +  diagSigmVec(k)*imresize(temp(cropM+1:end-cropM-1,cropN+1:end-cropN-1),[M N]);
end
signalFiltered2D  = signalRevRot2D/sum(diagSigmVec);