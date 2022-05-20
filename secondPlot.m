function [plot2] = secondPlot(fck,leq,b)
[sigmac2,ech2] = secondCompressionBlock(fck,leq,b);
[sigmat2,etck] = secondTensionBlock(fck,leq);

plot2=figure(2);
subplot(2,1,1)
plot(ech2,sigmac2,'linewidth',2)
xlabel('Deformação Inelastica de Compressão','FontSize',14)
ylabel('Tensão de Compressão [MPa]','FontSize',14)
title('Gráfico Tensão x Deformação Inelastica de Compressão do Concreto','FontSize',16)

subplot(2,1,2)
plot(etck,sigmat2,'linewidth',2)
xlabel('Deformação Inelastica de Tração','FontSize',14)
ylabel('Tensão de Tração[MPa]','FontSize',14)
title('Gráfico Tensão x Deformação Inelastica de Tração do Concreto','FontSize',16)
end