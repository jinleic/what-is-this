︠3ae58ff3-0c36-4197-8f69-0ff80629e8c0s︠
def H(x):
    if x <=0:
        return 0
    if x>=1:
        return 0
    else:
        return -(x*log(x,2)+(1-x)*log(1-x,2))
︡ff334c03-bb66-4052-a292-2611c80eadac︡{"done":true}
︠5a2e52be-16cb-4b4b-9c8c-b5bc0f187bc6s︠
def g(x,y):
    return ((2-2*y)*H((1-x)/2)- H((1-x)*(1-y)))
︡1c5b52ce-d7a0-4dd5-8795-05f3ed78bab9︡{"done":true}
︠9c69fa0c-8524-4384-8646-d8f35cf7fc45s︠

x, y = var('x y')
derivative(g(x,y),x)
derivative(derivative(g(x,y),x),x)

︡2a67bdf5-1029-4c75-b3c9-01945aa227d4︡{"stdout":"(y - 1)*(log(1/2*x + 1/2)/log(2) - log(-1/2*x + 1/2)/log(2)) + (y - 1)*log((x - 1)*(y - 1))/log(2) - (y - 1)*log(-(x - 1)*(y - 1) + 1)/log(2)\n"}︡{"stdout":"(y - 1)*(1/((x + 1)*log(2)) - 1/((x - 1)*log(2))) - (y - 1)^2/(((x - 1)*(y - 1) - 1)*log(2)) + (y - 1)/((x - 1)*log(2))\n"}︡{"done":true}
︠d42e73ca-afbd-4708-84f7-21a960556e4c︠
︡3f8b7517-8b03-42a8-82f8-3491899788ee︡
︠25479f43-469b-45ac-a45e-58154c1b0991s︠
plot(2*g(0.39,y)-1.044*g(0,y),(y,0.5,1))
plot(2*g(0.39,y)-1.044*g(0,y),(y,0.5,0.50001))
plot(2*g(0.39,y)-1.044*g(0,y),(y,0.9999,1))
2*g(0.39,1)-1.044*g(0,1)
︡6e6e25c5-1a19-4d7a-a15f-604bc84867e8︡{"file":{"filename":"/tmp/tmpop84i6_q/tmp_ahvg8b4b.svg","show":true,"text":null,"uuid":"2b7e7dcb-1d34-4504-861d-b1befc3be5f8"},"once":false}︡{"file":{"filename":"/tmp/tmpop84i6_q/tmp_slm_c5t9.svg","show":true,"text":null,"uuid":"cc3e0536-d83f-4426-8777-1f694707b94b"},"once":false}︡{"file":{"filename":"/tmp/tmpop84i6_q/tmp_icjn76ub.svg","show":true,"text":null,"uuid":"84150ecb-25a8-40be-85fb-f8d51412c2ea"},"once":false}︡{"stdout":"0.000000000000000\n"}︡{"done":true}
︠77ce7b0e-0164-4071-afce-7864bc715456︠









