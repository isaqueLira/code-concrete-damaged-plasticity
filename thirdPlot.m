function [plot3] = thirdPlot(fck,leq,b)
[ecu] = firstCompressionBlock(fck,leq,b);
[ech2,dc2] = secondCompressionBlock(fck,leq,b);
[etu] = eqInitialTension(fck,leq);
[etck,dt] = secondTensionBlock(fck,leq);

plot3=figure(3);
subplot(2,1,1)
plot(ech2,dc2,'linewidth',2)
xlabel('Deforma��o Inelastica de Compress�o','FontSize',14)
ylabel('Dano em Compress�o','FontSize',14)
title('Gr�fico Dano x Deforma��o Inelastica de Compress�o do Concreto','FontSize',16)
axis([0 ecu 0 1.1])

subplot(2,1,2)
plot(etck,dt,'linewidth',2)
xlabel('Deforma��o Inelastica de Tra��o','FontSize',14)
ylabel('Dano em Tra��o','FontSize',14)
title('Gr�fico Dano x Deforma��o Inelastica de Tra��o do Concreto','FontSize',16)
axis([0 etu 0 1.1])
end