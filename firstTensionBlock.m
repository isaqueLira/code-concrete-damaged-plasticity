function [et1,sigmat1] = firstTensionBlock(fck,leq)
[ftm,E0] = eqInitialCompression(fck,leq);
[etm,etu] = eqInitialTension(fck,leq);

et1=linspace(0,etu,100); % Deforma��o de Tra��o
sigmat1=zeros(length(et1),1)'; % Tens�o de Tra��o
for i=1:length(et1) % Loop para a elabora��o das curvas tens�oxdeforma��o e danoxdeforma��o
    if et1(i)<=etm % Trecho Linear
        sigmat1(i)=et1(i)*E0; % Alfarah et al. (2017)
    else % Trecho N�o-Linear
        sigmat1(i)=(((1+(3*((et1(i)-etm)/(etu-etm)))^3)*(exp(-6.93*((et1(i)-etm)/(etu-etm)))))-(((et1(i)-etm)/(etu-etm))*(1+(3^3))*(exp(-6.93))))*ftm; % Adaptada de Hordijk (1992)
    end
end
end