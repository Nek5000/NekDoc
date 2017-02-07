import patterns;
pair z0=(0,0);
real r=30;
int n=20;
real h=0.1;

label("$\Omega_f$",(0,10));
path c2=ellipse(z0,r+20,r+5);
draw(c2);

add("hatch",hatch(1mm));
path c1=ellipse(z0,r-8,r);
draw(c1);
filldraw(c2,pattern("hatch"));

label("$\Omega_S$",(10,40));
filldraw(c1,white);
//real par(real x,real a,real b,real c) {return -(a*x^2+b*x+c);}
     
//real sum= -3;
//real prod=2;

//guide fv;

//for(int i=0; i < n; ++i) {	
//      real t=i*h;		
//      pair z=(t,par(t,1,sum,prod));
//      fv=fv..z;	
//}	
 
//draw("$u$",(shift+50,par(shift+50,1,sum,prod)/60+shift)--(shift+70,par(shift+70,1,sum,prod)/50+shift),red,EndArrow);     
//draw(fv);

//label("$0$",(0,10));

