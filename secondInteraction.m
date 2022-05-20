function [yc,ecu,ec,sigmac,ech,ecpl,dc,b_n] = secondInteraction(fck,leq,b)
[fcm,E0,Gch,ecel,ecm] = eqInitialCompression(fck,leq);
[yc,gc3,ecu,ec] = eqFinalCompression(fck,leq,b);
[sigmac,ech,dc,ecpl,b1_mean] = firstInteraction(fck,leq,b);

while abs(b1_mean-b)>=0.00001
    b=b1_mean;
    yc=((pi^2)*fcm*ecm)/(2*(((Gch/leq)-(0.5*fcm)*((ecm*(1-b))+(b*(fcm/E0))))^2)); % Eq. 28 de Alfarah et al. (2017) 
    ecu=((yc*ecm)+2*ecm*(tan(gc3*sqrt(yc/(2*fcm*ecm))))*(sqrt(yc/(2*fcm*ecm))))/yc; % Alfarah et al. (2017)
    ec=linspace(ecel,ecu,100); % Deformação total de compressão
    
    [sigmac,ech,dc,ecpl,b1_mean] = firstInteraction(fck,leq,b); %Function First Interaction
end
b_n=b1_mean;
end