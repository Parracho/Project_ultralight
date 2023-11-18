(*
     ==============================
     *  CalcHEP  3.8.9 *
     ==============================
  process  h2(p1)+h2(p2)->m(p3)+M(p4)
*)

parameters={
 EE -> 3.13330000000*10^(-1)
,Mm -> 1.05000000000*10^(-1)
,SW -> 4.74000000000*10^(-1)
,MW -> 8.03850000000*10^(1)
,a -> 2.30000000000*10^(1)
,Mh -> 1.00000000000*10^(0)
,mh2 -> 1.00000000000*10^(0)
,HHH -> 1.00000000000*10^(0)
,hHH -> 1.00000000000*10^(0)
,wh -> 0.00000000000*10^(0)
,wh2 -> 0.00000000000*10^(0)
           };

substitutions={
 cosa->Cos[a]
,sina->Sin[a]
              };

inParticles = {"h2","h2"}
outParticles = {"m","M"}

SetAttributes[ SC, Orderless ];

SC[ a_ , b_ + c_ ] := SC[a,b]+SC[a,c];

SC[ x_?NumberQ * a_ , b_ ] := x * SC[ a, b ]



p4 = +p1+p2-p3;
p1/: SC[p1,p1] =mh2^2;
p2/: SC[p2,p2] =mh2^2;
p3/: SC[p3,p3] =Mm^2;
p2/: SC[p2,p3] = -1*(Mm^2-mh2^2-mh2^2-Mm^2-2*SC[p1,p2]+2*SC[p1,p3])/2;

initSum;

(*
  Diagram  1 in subprocess 10
*)
totFactor = ((cosa^2*hHH^2*Mm^2*EE^2)/(MW^2*SW^2));
numerator =(SC[p1,p2]+mh2^2-2*Mm^2);
denominator =(propDen[-p1-p2,Mh,wh]^2);

addToSum;

(*
  Diagram  2 in subprocess 10
*)
totFactor = ((-2*cosa*sina*hHH*HHH*Mm^2*EE^2)/(MW^2*SW^2));
numerator =(SC[p1,p2]+mh2^2-2*Mm^2);
denominator =(propDen[-p1-p2,Mh,wh]*propDen[-p1-p2,mh2,wh2]);

addToSum;

(*
  Diagram  3 in subprocess 10
*)
totFactor = ((sina^2*HHH^2*Mm^2*EE^2)/(MW^2*SW^2));
numerator =(SC[p1,p2]+mh2^2-2*Mm^2);
denominator =(propDen[-p1-p2,mh2,wh2]^2);

addToSum;

finishSum;
