//pair z0=(0,0);
//real r=30;
//path c1=circle(z0,r);
//draw(c1);

//draw((-100,-100)--(-100,50)); 
//for (int i=1; i<31; ++i)
//{
//  draw((-105,-100+(i-1)*5)--(-100,-100+i*5)); 
//}

//draw((100,-100)--(100,50)); 
//for (int i=1; i<31; ++i)
//{
//  draw((100,-100+(i)*5)--(105,-100+(i-1)*5)); 
//}
//draw("$H$",(0,0)--(100,0),darkgreen,Arrows,Bars,PenMargins);
//draw("$\ell=2H$",(-100,-60)--(100,-60),darkgreen,Arrows,Bars,PenMargins);
//draw("$r$",(0,0)--(14,14));
//draw("$g$",(-40,0)--(-40,-40),EndArrow);
//size(200);
//real par1(real x,real b,real c) {return -(x^2+b*x+c);}
real par(real x,real a,real b,real c) {return -(a*x^2+b*x+c);}
int n=20;
real a=1.5;
real width=2;
real shift=150; 
     
real sum= -(2*shift+100);
real prod=shift*(shift+100);

guide fv;
guide gu;
     
for(int i=0; i < n; ++i) {
      real t=shift+i*width;
      
      pair z=(t,par(t,1,sum,prod)/60+shift);
//      draw((t,par(t,1,sum,prod)/60+shift)--(t+10,par(t+10,2,sum,0)/60+shift),EndArrow); 
      fv=fv..z;
      pair y=(t+38,par(t+38,1,sum,prod)/60+shift);
      gu=gu..y;
}
 
draw("$v$",(shift+10,par(shift+10,1,sum,prod)/60+shift)--(shift+30,par(shift+30,1,sum,prod)/50+shift),blue,EndArrow);     
draw("$u$",(shift+50,par(shift+50,1,sum,prod)/60+shift)--(shift+70,par(shift+70,1,sum,prod)/50+shift),red,EndArrow);     
draw(fv);

draw(gu,red);
draw((0,0)--(0,200),EndArrow);
draw((0,0)--(200,0),EndArrow);

draw((150,150)--(150,200),EndArrow);
draw((150,150)--(200,150),EndArrow);
label("$0$",(shift,shift-10));
label("$v(0)=u(t)$",(shift,shift-20));
label("$s$",(shift+38,shift-10));
real t=shift;
label("$Q_{f}^{t^*}u$",(t+38,par(t+38,1,sum,prod)/60+shift+10));
label("$*$",(t+38,par(t+38,1,sum,prod)/60+shift));
label("$t$",(shift,-10));
label("$t^*$",(shift+38,-10));
label("$*$",(shift,0));
label("$*$",(shift+38,0));
