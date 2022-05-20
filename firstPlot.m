function [plot1] = firstPlot(fck,leq,b)
[ec1,sigmac1] = firstCompressionBlock(fck,leq,b);
[et1,sigmat1] = firstTensionBlock(fck,leq);

plot1=figure(1);
subplot(2,1,1)
plot(ec1,sigmac1,'linewidth',2)
xlabel('Deforma��o de Compress�o','FontSize',14)
ylabel('Tens�o de Compress�o [MPa]','FontSize',14)
title('Gr�fico Tens�o x Deforma��o de Compress�o do Concreto','FontSize',16)

subplot(2,1,2)
plot(et1,sigmat1,'linewidth',2)
xlabel('Deforma��o de Tra��o','FontSize',14)
ylabel('Tens�o de Tra��o[MPa]','FontSize',14)
title('Gr�fico Tens�o x Deforma��o de Tra��o do Concreto','FontSize',16)
end