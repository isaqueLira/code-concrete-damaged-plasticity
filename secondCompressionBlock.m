function [ec2,sigmac2,ech2,ecpl2,dc2,b2] = secondCompressionBlock(fck,leq,b)
[fcm,E0,ecm] = eqInitialCompression(fck,leq);
[ac,bc] = eqFinalCompression(fck,leq,b);
[yc,ecu] = firstCompressionBlock(fck,leq,b);

ec2=linspace(ecm,ecu,50); % Deforma��o total de compress�o do trecho n�o-linear
sigmac2=zeros(length(ec2),1)'; % Tens�o de compress�o do trecho n�o-linear
ech2=zeros(length(ec2),1)'; % Deforma��o inelastica de compress�o
ecpl2=zeros(length(ec2),1)'; % Deforma��o pl�stica de compress�o
dc2=zeros(length(ec2),1)'; % Dano de compress�o
b2=zeros(length(ec2),1)'; % Raz�o deforma��o pl�stica de compress�o e deforma��o inelastica de compress�o 
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