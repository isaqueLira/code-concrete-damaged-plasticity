function [plot2] = secondPlot(fck,leq,b)
[sigmac2,ech2] = secondCompressionBlock(fck,leq,b);
[sigmat2,etck] = secondTensionBlock(fck,leq);

plot2=figure(2);
subplot(2,1,1)
plot(ech2,sigmac2,'linewidth',2)
xlabel('Deforma��o Inelastica de Compress�o','FontSize',14)
ylabel('Tens�o de Compress�o [MPa]','FontSize',14)
title('Gr�fico Tens�o x Deforma��o Inelastica de Compress�o do Concreto','FontSize',16)

subplot(2,1,2)
plot(etck,sigmat2,'linewidth',2)
xlabel('Deforma��o Inelastica de Tra��o','FontSize',14)
ylabel('Tens�o de Tra��o[MPa]','FontSize',14)
title('Gr�fico Tens�o x Deforma��o Inelastica de Tra��o do Concreto','FontSize',16)
end