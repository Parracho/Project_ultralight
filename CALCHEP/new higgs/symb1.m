(*
     ==============================
     *  CalcHEP  3.8.9 *
     ==============================
  process  h2(p1)+h2(p2)->~dm(p3)+~dm(p4)
*)

parameters={
 vsigma -> 1.00000000000*10^(0)
,mius -> 1.00000000000*10^(0)
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
,Mdm->Sqrt[2*mius^2]
              };

inParticles = {"h2","h2"}
outParticles = {"~dm","~dm"}

SetAttributes[ SC, Orderless ];

SC[ a_ , b_ + c_ ] := SC[a,b]+SC[a,c];

SC[ x_?NumberQ * a_ , b_ ] := x * SC[ a, b ]



p4 = +p1+p2-p3;
p1/: SC[p1,p1] =mh2^2;
p2/: SC[p2,p2] =mh2^2;
p3/: SC[p3,p3] =Mdm^2;
p2/: SC[p2,p3] = -1*(Mdm^2-mh2^2-mh2^2-Mdm^2-2*SC[p1,p2]+2*SC[p1,p3])/2;

initSum;

(*
  Diagram  1 in subprocess 1
*)
totFactor = ((sina^2*hHH^2*Mh^4)/(2*vsigma^2));
numerator =(1);
denominator =(propDen[-p1-p2,Mh,wh]^2);

addToSum;

(*
  Diagram  2 in subprocess 1
*)
totFactor = ((cosa*sina*hHH*HHH*mh2^2*Mh^2)/(vsigma^2));
numerator =(1);
denominator =(propDen[-p1-p2,Mh,wh]*propDen[-p1-p2,mh2,wh2]);

addToSum;

(*
  Diagram  3 in subprocess 1
*)
totFactor = ((cosa^2*HHH^2*mh2^4)/(2*vsigma^2));
numerator =(1);
denominator =(propDen[-p1-p2,mh2,wh2]^2);

addToSum;

finishSum;
