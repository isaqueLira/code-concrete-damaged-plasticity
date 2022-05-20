function [plot1] = firstPlot(fck,leq,b)
[ec1,sigmac1] = firstCompressionBlock(fck,leq,b);
[et1,sigmat1] = firstTensionBlock(fck,leq);

plot1=figure(1);
subplot(2,1,1)
plot(ec1,sigmac1,'linewidth',2)
xlabel('Deformação de Compressão','FontSize',14)
ylabel('Tensão de Compressão [MPa]','FontSize',14)
title('Gráfico Tensão x Deformação de Compressão do Concreto','FontSize',16)

subplot(2,1,2)
plot(et1,sigmat1,'linewidth',2)
xlabel('Deformação de Tração','FontSize',14)
ylabel('Tensão de Tração[MPa]','FontSize',14)
title('Gráfico Tensão x Deformação de Tração do Concreto','FontSize',16)
end