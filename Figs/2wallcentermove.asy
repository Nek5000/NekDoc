pair z0=(0,0);
real r=30;
path c1=circle(z0,r);
draw(c1);

draw((-100,-100)--(-100,50)); 
for (int i=1; i<31; ++i)
{
  draw((-105,-100+(i-1)*5)--(-100,-100+i*5)); 
}

draw((100,-100)--(100,50)); 
for (int i=1; i<31; ++i)
{
  draw((100,-100+(i)*5)--(105,-100+(i-1)*5)); 
}
draw("$H$",(0,0)--(100,0),darkgreen,Arrows,Bars,PenMargins);
draw("$\ell=2H$",(-100,-60)--(100,-60),darkgreen,Arrows,Bars,PenMargins);
draw("$r$",(0,0)--(14,14));
draw("$g$",(-40,0)--(-40,-40),EndArrow);


draw((-80,70)--(-140,70),EndArrow);
draw((80,70)--(140,70),EndArrow);
