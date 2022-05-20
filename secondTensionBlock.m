function [et2,sigmat2,etck,etpl,dt,b3] = secondTensionBlock(fck,leq)
[ftm,E0] = eqInitialCompression(fck,leq);
[etm,etu,at,bt] = eqInitialTension(fck,leq);

et2=linspace(etm,etu,50); % Deforma��o de Tra��o
sigmat2=zeros(length(et2),1)'; % Tens�o de Tra��o
etck=zeros(length(et2),1)'; % Deforma��o Inelastica de Tra��o
etpl=zeros(length(et2),1)'; % Deforma��o Pl�stica de Tra��o
dt=zeros(length(et2),1)'; % Dano de Tra��o
b3=zeros(length(et2),1)';
for i=1:length(et2) % Loop para a elabora��o das curvas tens�oxdeforma��o e danoxdeforma��o
    sigmat2(i)=(((1+(3*((et2(i)-etm)/(etu-etm)))^3)*(exp(-6.93*((et2(i)-etm)/(etu-etm)))))-(((et2(i)-etm)/(etu-etm))*(1+(3^3))*(exp(-6.93))))*ftm; % Adaptada de Hordijk (1992)
    etck(i)=et2(i)-(sigmat2(i)/E0); % Adaptada de Birtel e Mark (2006) e Alfarah et al. (2017)
    dt(i)=1-((1/(2+at))*(((2*(1+at))*exp(-bt*etck(i)))-(at*exp((-2*bt)*etck(i)))));
    etpl(i)=etck(i)-((sigmat2(i)/E0)*(dt(i)/(1-dt(i))));
    b3(i)=etpl(i)/etck(i);
end
end