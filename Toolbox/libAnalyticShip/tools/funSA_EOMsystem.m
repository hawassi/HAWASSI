function FX=funSA_EOMsystem(XX,ww,M,I,addMass_mat,dampCoef_mat,Cvect,F0ExcVect,PhaseExcVect,bnu2,bnu4)
a22=addMass_mat(1,1);a23=addMass_mat(1,2);a24=addMass_mat(1,3);
a32=addMass_mat(2,1);a33=addMass_mat(2,2);a34=addMass_mat(2,3);
a42=addMass_mat(3,1);a43=addMass_mat(3,2);a44=addMass_mat(3,3);

b22=dampCoef_mat(1,1);b23=dampCoef_mat(1,2);b24=dampCoef_mat(1,3);
b32=dampCoef_mat(2,1);b33=dampCoef_mat(2,2);b34=dampCoef_mat(2,3);
b42=dampCoef_mat(3,1);b43=dampCoef_mat(3,2);b44=dampCoef_mat(3,3);

c22=Cvect(1);
c33=Cvect(2);
c44=Cvect(3);
F02=F0ExcVect(1);
F03=F0ExcVect(2);
F04=F0ExcVect(3);
varphif2=PhaseExcVect(1);
varphif3=PhaseExcVect(2);
varphif4=PhaseExcVect(3);


X2=XX(1);X3=XX(2);X4=XX(3);
varphix2=XX(4);varphix3=XX(5);varphix4=XX(6);

varphifx2=varphif2+varphix2;
varphifx3=varphif3+varphix3;
varphifx4=varphif4+varphix4;

FX(1)=(M+a22)*ww^2*sin(varphix2)*X2+(b22+bnu2*abs(ww*cos(varphix2)*X2))*ww*cos(varphix2)*X2-c22*sin(varphix2)*X2-a23*ww^2*sin(varphif2-varphifx3)*X3+b23*ww*cos(varphif2-varphifx3)*X3-a24*ww^2*sin(varphif2-varphifx4)*X4+b24*ww*cos(varphif2-varphifx4)*X4;
FX(2)=(M+a33)*ww^2*sin(varphix3)*X3+b33*ww*cos(varphix3)*X3-c33*sin(varphix3)*X3-a34*ww^2*sin(varphif3-varphifx4)*X4+b34*ww*cos(varphif3-varphifx4)*X4-a32*ww^2*sin(varphif3-varphifx2)*X2+b32*ww*cos(varphif3-varphifx2)*X2;
FX(3)=(I+a44)*ww^2*sin(varphix4)*X4+(b44+bnu4*abs(ww*cos(varphix4)*X4))*ww*cos(varphix4)*X4-c44*sin(varphix4)*X4-a43*ww^2*sin(varphif4-varphifx3)*X3+b43*ww*cos(varphif4-varphifx3)*X3-a42*ww^2*sin(varphif4-varphifx2)*X2+b42*ww*cos(varphif4-varphifx2)*X2;%
FX(4)=-(M+a22)*ww^2*cos(varphix2)*X2+(b22+bnu2*abs(ww*sin(varphix2)*X2))*ww*sin(varphix2)*X2+c22*cos(varphix2)*X2-F02 -a23*ww^2*cos(varphif2-varphifx3)*X3-b23*ww*sin(varphif2-varphifx3)*X3-a24*ww^2*cos(varphif2-varphifx4)*X4-b24*ww*sin(varphif2-varphifx4)*X4;
FX(5)=-(M+a33)*ww^2*cos(varphix3)*X3+b33*ww*sin(varphix3)*X3+c33*cos(varphix3)*X3-F03-a34*ww^2*cos(varphif3-varphifx4)*X4-b34*ww*sin(varphif3-varphifx4)*X4-a32*ww^2*cos(varphif3-varphifx2)*X2-b32*ww*sin(varphif3-varphifx2)*X2;
FX(6)=-(I+a44)*ww^2*cos(varphix4)*X4+(b44+bnu4*abs(ww*sin(varphix4)*X4))*ww*sin(varphix4)*X4+c44*cos(varphix4)*X4-F04-a43*ww^2*cos(varphif4-varphifx3)*X3-b43*ww*sin(varphif4-varphifx3)*X3-a42*ww^2*cos(varphif4-varphifx2)*X2-b42*ww*sin(varphif4-varphifx2)*X2;%