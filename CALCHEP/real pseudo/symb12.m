(*
     ==============================
     *  CalcHEP  3.8.9 *
     ==============================
  process  ~dm(p1)+~dm(p2)->Z(p3)+Z(p4)
*)

parameters={
 EE -> 3.13330000000*10^(-1)
,SW -> 4.74000000000*10^(-1)
,MW -> 8.03850000000*10^(1)
,vsigma -> 1.00000000000*10^(0)
,mius -> 1.00000000000*10^(0)
,a -> 2.30000000000*10^(1)
,Mh -> 1.00000000000*10^(0)
,mh2 -> 1.00000000000*10^(0)
,wh -> 0.00000000000*10^(0)
,wh2 -> 0.00000000000*10^(0)
           };

substitutions={
 cosa->Cos[a]
,sina->Sin[a]
,Mdm->Sqrt[2*mius^2]
,MZ->MW/CW
,CW->Sqrt[1-SW^2]
              };

inParticles = {"~dm","~dm"}
outParticles = {"Z","Z"}

SetAttributes[ SC, Orderless ];

SC[ a_ , b_ + c_ ] := SC[a,b]+SC[a,c];

SC[ x_?NumberQ * a_ , b_ ] := x * SC[ a, b ]



p4 = +p1+p2-p3;
p1/: SC[p1,p1] =Mdm^2;
p2/: SC[p2,p2] =Mdm^2;
p3/: SC[p3,p3] =MZ^2;
p2/: SC[p2,p3] = -1*(MZ^2-Mdm^2-Mdm^2-MZ^2-2*SC[p1,p2]+2*SC[p1,p3])/2;

initSum;

(*
  Diagram  1 in subprocess 12
*)
totFactor = ((cosa^2*sina^2*mh2^4*MW^2*EE^2)/(2*MZ^4*CW^4*vsigma^2*SW^2));
numerator =(SC[p1,p2]^2+2*SC[p1,p2]*Mdm^2-2*SC[p1,p2]*MZ^2+Mdm^4-2*Mdm^2*MZ^
 2+3*MZ^4);
denominator =(propDen[-p1-p2,mh2,wh2]^2);

addToSum;

(*
  Diagram  2 in subprocess 12
*)
totFactor = ((-cosa^2*sina^2*mh2^2*Mh^2*MW^2*EE^2)/(MZ^4*CW^4*vsigma^2*SW^2
 ));
numerator =(SC[p1,p2]^2+2*SC[p1,p2]*Mdm^2-2*SC[p1,p2]*MZ^2+Mdm^4-2*Mdm^2*MZ^
 2+3*MZ^4);
denominator =(propDen[-p1-p2,mh2,wh2]*propDen[-p1-p2,Mh,wh]);

addToSum;

(*
  Diagram  3 in subprocess 12
*)
totFactor = ((cosa^2*sina^2*Mh^4*MW^2*EE^2)/(2*MZ^4*CW^4*vsigma^2*SW^2));
numerator =(SC[p1,p2]^2+2*SC[p1,p2]*Mdm^2-2*SC[p1,p2]*MZ^2+Mdm^4-2*Mdm^2*MZ^
 2+3*MZ^4);
denominator =(propDen[-p1-p2,Mh,wh]^2);

addToSum;

finishSum;
