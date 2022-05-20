function [yc,gc2,gc3,ecu,ec,ac,bc] = eqFinalCompression(fck,leq,b)
[fcm,fc0,E0,Gch,gc,ecel,ecm] = eqInitialCompression(fck,leq);

yc=((pi^2)*fcm*ecm)/(2*(((Gch/leq)-(0.5*fcm)*((ecm*(1-b))+(b*(fcm/E0))))^2)); % Eq. 28 de Alfarah et al. (2017) 
gc2=((ecm-ecel)*fcm)/1.2; % Krätzig e Pölling (2004) 
gc3=gc-gc2; % Krätzig e Pölling (2004)
ecu=((yc*ecm)+2*ecm*(tan(gc3*sqrt(yc/(2*fcm*ecm))))*(sqrt(yc/(2*fcm*ecm))))/yc; % Alfarah et al. (2017)
ec=linspace(ecel,ecu,300); % Deformação total de compressão
ac=((2*(fcm/fc0))-1)+(2*sqrt(((fcm/fc0)^2)-(fcm/fc0))); % Eq. 24 de Alfarah et al. (2017)
bc=((1.97*(fck+8))/Gch)*leq; % Eq. 18 de Alfarah et. al (2017)
end