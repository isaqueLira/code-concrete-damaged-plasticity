function [ec2,sigmac2,ech2,ecpl2,dc2,b2] = secondCompressionBlock(fck,leq,b)
[fcm,E0,ecm] = eqInitialCompression(fck,leq);
[ac,bc] = eqFinalCompression(fck,leq,b);
[yc,ecu] = firstCompressionBlock(fck,leq,b);

ec2=linspace(ecm,ecu,50); % Deformação total de compressão do trecho não-linear
sigmac2=zeros(length(ec2),1)'; % Tensão de compressão do trecho não-linear
ech2=zeros(length(ec2),1)'; % Deformação inelastica de compressão
ecpl2=zeros(length(ec2),1)'; % Deformação plástica de compressão
dc2=zeros(length(ec2),1)'; % Dano de compressão
b2=zeros(length(ec2),1)'; % Razão deformação plástica de compressão e deformação inelastica de compressão 
for i=1:length(ec2)
    if ec2(i)==ecm
        sigmac2(i)=(((2+(yc*fcm*ecm))/(2*fcm))-(ec2(i)*yc)+((ec2(i)^2)*(yc/(2*ecm))))^-1;
        ech2(i)=0;
        dc2(i)=0;
        ecpl2(i)=0;
    else
        sigmac2(i)=(((2+(yc*fcm*ecm))/(2*fcm))-(ec2(i)*yc)+((ec2(i)^2)*(yc/(2*ecm))))^-1;
        ech2(i)=ec2(i)-(sigmac2(i)/E0);
        dc2(i)=1-((1/(2+ac))*(((2*(1+ac))*exp(-bc*ech2(i)))-(ac*exp((-2*bc)*ech2(i))))); % Eq. 19 de Alfarah et. al (2017)
        if dc2(i)<=0.994
            ecpl2(i)=ech2(i)-((sigmac2(i)/E0)*(dc2(i)/(1-dc2(i))));
            b2(i)=ecpl2(i)/ech2(i);
        else
            dc2(i)=0.994;
            ecpl2(i)=ech2(i)-((sigmac2(i)/E0)*(dc2(i)/(1-dc2(i))));
            b2(i)=ecpl2(i)/ech2(i);
        end
    end
end
end