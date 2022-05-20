function [b,yc,ecu,ec1,sigmac1] = firstCompressionBlock(fck,leq,b)
[fcm,Eci,E0,Gch,ecel,ecm] = eqInitialCompression(fck,leq);
[gc3] = eqFinalCompression(fck,leq,b);
[b_n] = secondInteraction(fck,leq,b);

b=b_n;
yc=((pi^2)*fcm*ecm)/(2*(((Gch/leq)-(0.5*fcm)*((ecm*(1-b))+(b*(fcm/E0))))^2)); % Eq. 28 de Alfarah et al. (2017) 
ecu=((yc*ecm)+2*ecm*(tan(gc3*sqrt(yc/(2*fcm*ecm))))*(sqrt(yc/(2*fcm*ecm))))/yc; % Alfarah et al. (2017)
ec1=linspace(0,ecu,300); % Deformação total de compressão
sigmac1=zeros(length(ec1),1)'; % Tensão de compressão

for i=1:length(ec1)
    if ec1(i)<=ecel
        sigmac1(i)=round(ec1(i)*E0,2);
    elseif ec1(i)>ecel && ec1(i)<=ecm
        sigmac1(i)=round((((Eci*(ec1(i)/fcm)-((ec1(i)/ecm)^2)))/(1+(((Eci*(ecm/fcm))-2)*(ec1(i)/ecm))))*fcm,2);
    else
        sigmac1(i)=round((((2+(yc*fcm*ecm))/(2*fcm))-(ec1(i)*yc)+((ec1(i)^2)*(yc/(2*ecm))))^-1,2);
    end
end
end
