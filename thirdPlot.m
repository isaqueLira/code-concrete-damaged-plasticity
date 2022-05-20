function [plot3] = thirdPlot(fck,leq,b)
[ecu] = firstCompressionBlock(fck,leq,b);
[ech2,dc2] = secondCompressionBlock(fck,leq,b);
[etu] = eqInitialTension(fck,leq);
[etck,dt] = secondTensionBlock(fck,leq);

plot3=figure(3);
subplot(2,1,1)
plot(ech2,dc2,'linewidth',2)
xlabel('Deformação Inelastica de Compressão','FontSize',14)
ylabel('Dano em Compressão','FontSize',14)
title('Gráfico Dano x Deformação Inelastica de Compressão do Concreto','FontSize',16)
axis([0 ecu 0 1.1])

subplot(2,1,2)
plot(etck,dt,'linewidth',2)
xlabel('Deformação Inelastica de Tração','FontSize',14)
ylabel('Dano em Tração','FontSize',14)
title('Gráfico Dano x Deformação Inelastica de Tração do Concreto','FontSize',16)
axis([0 etu 0 1.1])
end