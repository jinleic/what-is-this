︠3ae58ff3-0c36-4197-8f69-0ff80629e8c0︠
#first we define the function H, which is the binary entropy function
def H(x):
    if x <=0:
        return 0
    if x>=1:
        return 0
    else:
        return -(x*log(x,2)+(1-x)*log(1-x,2))
︡8ce36c56-64de-4dd6-b904-62ea561befe9︡{"done":true}
︠9fbfbebc-de08-4b20-93d9-8e09915b1db4︠
#we determine the values for a,b,c which we claim to be optimal
def ff(x):
    return(H(x)*(2-H(x))-H(2*x-x^2))
a=0.3
b=0.4
for n in range(0,1000):
    if ff((a+b)/2)<0:
        a=round((a+b)/2,100)
    else: 
        b=round((a+b)/2,100)
round(ff(b),200)
b
round(b,200)
a=(1-H(b))/(2-H(b))
c=a+(1-a)*b
round(c,100)
plot(ff,(0,1))
︡70cb2a19-ad4c-4923-a8c3-dab771424765︡{"stdout":"0.0\n"}︡{"stdout":"0.3294547385030371\n"}︡{"stdout":"0.3294547385030371\n"}︡{"stdout":"0.38234553336670285\n"}︡{"file":{"filename":"/tmp/tmpwcvkjkai/tmp_ut627shr.svg","show":true,"text":null,"uuid":"92148834-9d00-4650-9c92-972c0162d6bb"},"once":false}︡{"done":true}
︠8173f28d-ac62-4ac8-8b4e-dfe43173ed22︠
%we compute derivatives and set them as functions to compute the optimal choice for alpha

derivative((1-(c-x)/(1-x))^2*H(2*x-x^2)-(1-(c-x)/(1-x))*H(x),x)
︡4128eba7-88ca-414a-909b-133af3d8d55b︡{"stdout":"-2*((x - 1)*log(x^2 - 2*x + 1)/log(2) - (x - 1)*log(-x^2 + 2*x)/log(2))*((x - 0.38234553336670285)/(x - 1) - 1)^2 + 2*((x^2 - 2*x + 1)*log(x^2 - 2*x + 1)/log(2) - (x^2 - 2*x)*log(-x^2 + 2*x)/log(2))*((x - 0.38234553336670285)/(x - 1) - 1)*((x - 0.38234553336670285)/(x - 1)^2 - 1/(x - 1)) - ((x - 0.38234553336670285)/(x - 1) - 1)*(log(x)/log(2) - log(-x + 1)/log(2)) + (x*log(x)/log(2) - (x - 1)*log(-x + 1)/log(2))*((x - 0.38234553336670285)/(x - 1)^2 - 1/(x - 1))\n"}︡{"done":true}
︠35970402-67e6-4c31-905e-26cd1e47c411s︠
def g1(x):
    return( -2*((x - 1)*log(x^2 - 2*x + 1)/log(2) - (x - 1)*log(-x^2 + 2*x)/log(2))*((x - c)/(x - 1) - 1)^2 + 2*((x^2 - 2*x + 1)*log(x^2 - 2*x + 1)/log(2) - (x^2 - 2*x)*log(-x^2 + 2*x)/log(2))*((x - c)/(x - 1) - 1)*((x - c)/(x - 1)^2 - 1/(x - 1)) - ((x - c)/(x - 1) - 1)*(log(x)/log(2) - log(-x + 1)/log(2)) + (x*log(x)/log(2) - (x - 1)*log(-x + 1)/log(2))*((x - c)/(x - 1)^2 - 1/(x - 1)))
︡f28d35b6-8409-4b14-b011-70f276e48995︡{"done":true}
︠3d41a264-1aba-4509-8f6d-4087958adfee︠
derivative((1-2*(c-x)/(1-x))-(1-(c-x)/(1-x))*H(x),x)
︡4bf4f625-c887-408e-9f15-694bb5f1795a︡{"stdout":"-((x - 0.38234553336670285)/(x - 1) - 1)*(log(x)/log(2) - log(-x + 1)/log(2)) + (x*log(x)/log(2) - (x - 1)*log(-x + 1)/log(2))*((x - 0.38234553336670285)/(x - 1)^2 - 1/(x - 1)) + (2*x - 0.7646910667334057)/(x - 1)^2 - 2/(x - 1)\n"}︡{"done":true}
︠44eee256-200e-41d8-8032-a502313ffd6fs︠
def g2(x):
    return(-((x - c)/(x - 1) - 1)*(log(x)/log(2) - log(-x + 1)/log(2)) + (x*log(x)/log(2) - (x - 1)*log(-x + 1)/log(2))*((x - c)/(x - 1)^2 - 1/(x - 1)) + (2*x - 2*c)/(x - 1)^2 - 2/(x - 1))
︡2daf7537-c9a2-4e12-b794-e3e3aa9e57c8︡{"done":true}
︠bb093565-b756-420e-bc9e-279858d52ea0s︠

alpha=-g1(b)/(-g1(b)+g2(b))
round(alpha,30)

︡db0b40b1-046a-4096-87b1-49671a7ee31f︡{"stdout":"0.035606980437447984\n"}︡{"done":true}
︠94e6358e-b83f-4d98-9641-ccf8a1b987f4s︠
#we first verify that the minimum of the parabola for fixed b is always larger than the maximum allowed value for a
def g(b):
    return (((1-alpha)*2*H(2*b-b^2)+alpha*2*H(min(2*b,0.5))-H(b))/((1-alpha)*2*H(2*b-b^2)) -(c-b)/(1-b))
plot(g,(0,c))
︡d7fdd147-d83a-4c0e-9389-69e83d0dd9f5︡{"file":{"filename":"/tmp/tmpwcvkjkai/tmp_uhlr58k9.svg","show":true,"text":null,"uuid":"185b7dc4-40b5-4c90-96d6-bb5a7d7abbfe"},"once":false}︡{"done":true}
︠aaa18f37-bdea-43aa-8f2a-8c0d0a62ccb0s︠
#finally, we verify that the function under consideration is positive at its minimum, i.e., for a=(c-b)/(1-b)
def f(a,b):
    return((1-alpha)*(1-a)^2*H(2*b-b^2)+alpha*(1-2*a)*H(min(2*b,0.5))-(1-a)*H(b))
def g(b):
    return f((c-b)/(1-b),b)

plot(g,(0.,c))
plot(g,(0.329,0.33))

︡d1bde7b9-8111-4b60-85cd-8cf3f2e4d01c︡{"file":{"filename":"/tmp/tmpwcvkjkai/tmp_ue4pkr9f.svg","show":true,"text":null,"uuid":"5e32149c-6d31-40e4-8fc6-e6cb36b1f90c"},"once":false}︡{"file":{"filename":"/tmp/tmpwcvkjkai/tmp_0i2jxlbu.svg","show":true,"text":null,"uuid":"3526935b-547c-444c-98cd-c03d22be024c"},"once":false}︡{"done":true}









