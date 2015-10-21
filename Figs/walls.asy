
import patterns;
pair z0=(70,0);
real r=30;
int n=20;
real h=0.1;
real x=-100;
real dx=20;

draw((-100,-100)--(200,-100)); 
for (int i=1; i<61; ++i)
{
  draw((-100+(i-1)*5,-105)--(-100+i*5,-100)); 
}
add("brick",brick(3mm,grey));
path r1=((x,50+0) -- (x,50+dx) -- (x+300,50+dx+10) -- (x+300,50+0) -- cycle);
path r2=((-100,-100) -- (-100,50) -- (200,50) -- (200,-100) -- cycle);
draw(r1);
filldraw(r1,pattern("hatch"));
filldraw(r2,paleblue);


path c2=ellipse(z0,r+15,r+2);
draw(c2);

add("hatch",hatch(1mm,grey));
path c1=ellipse(z0,r-8,r);
//draw(c1);
filldraw(c2,white);
filldraw(c2,pattern("hatch"));
//filldraw(c1,paleblue);

label("$\mathbf{\Omega_S}$",(80,60));
label("$\mathbf{\Omega_S}$",(100,0));
//label("$\mathbf{\Omega_f}$",(70,0));
label("$\mathbf{\Omega_f}$",(150,0));

draw("$\mathbf{\partial\overline{\Omega_f}(BC_1)}$",(-100,-40)--(-100,0),darkgreen,Arrows,PenMargins);
draw("$\mathbf{\partial\overline{\Omega_f}(BC_2)}$",(200,-40)--(200,0),darkgreen,Arrows,PenMargins);
draw((0,10)--(30,-100),darkgreen,EndArrow);
draw((0,10)--(20,50),darkgreen,EndArrow);
//draw((0,10)--(50,10),darkgreen,EndArrow);
draw((0,10)--(20,0),darkgreen,EndArrow);
label("$\mathbf{\partial\Omega_f=\partial\Omega_S}$",(-30,10));
draw("$\mathbf{\partial\overline{\Omega_S}}$",(60,100)--(60,80),EndArrow);

