function [sigmac,ech,dc,ecpl,b1,b1_mean] = firstInteraction(fck,leq,b)
[fcm,Eci,E0,ecel,ecm] = eqInitialCompression(fck,leq);
[yc,ec,ac,bc] = eqFinalCompression(fck,leq,b);

sigmac=zeros(length(ec),1)'; % Tensão de compressão
ech=zeros(length(ec),1)'; % Deformação inelastica de compressão
ecpl=zeros(length(ec),1)'; % Deformação plástica de compressão
dc=zeros(length(ec),1)'; % Dano de compressão
b1=zeros(length(ec),1)'; % Razão deformação plástica de compressão e deformação inelastica de compressão

for i=1:length(ec)
    if (ec(i)>=ecel) && (ec(i)<=ecm)
           sigmac(i)=round((((Eci*(ec(i)/fcm)-((ec(i)/ecm)^2)))/(1+(((Eci*(ecm/fcm))-2)*(ec(i)/ecm))))*fcm,3); % Tensão de compressão, fase 2
           ech(i)=ec(i)-(sigmac(i)/E0); 
           dc(i)=round(1-((1/(2+ac))*(((2*(1+ac))*exp(-bc*ech(i)))-(ac*exp((-2*bc)*ech(i))))),3); % Eq. 19 de Alfarah et. al (2017)
    elseif ec(i)>=ecm
        sigmac(i)=round((((2+(yc*fcm*ecm))/(2*fcm))-(ec(i)*yc)+((ec(i)^2)*(yc/(2*ecm))))^-1,3); % Tensão de compressão, fase 3
        ech(i)=ec(i)-(sigmac(i)/E0);
        dc(i)=round(1-((1/(2+ac))*(((2*(1+ac))*exp(-bc*ech(i)))-(ac*exp((-2*bc)*ech(i))))),3); % Eq. 19 de Alfarah et. al (2017)
    end
    if dc(i)<=0.994 % Valor estimado para garantir deformações plasticas sempre positivas
        ecpl(i)=ech(i)-((sigmac(i)/E0)*(dc(i)/(1-dc(i))));
    else
        dc(i)=0.994; % Valor estimado para garantir deformações plasticas sempre positivas
        ecpl(i)=ech(i)-((sigmac(i)/E0)*(dc(i)/(1-dc(i))));
    end
    b1(i)=ecpl(i)/ech(i);
end

b1_mean=mean(b1,'omitnan');
end