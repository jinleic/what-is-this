%https://www.mathworks.com/help/optim/ug/solve-nonlinear-optimization-problem-based.html
%https://www.mathworks.com/help/optim/ug/optimproblem.html
clc;
%beta=0.5;
%c=0.39;
%c=0.38196601125; accuaracy about 10^{-7}
%c=0.3823455

%C= linspace(0.382,0.383,5); 
%C=[0.3825,0.3828];% next: 0.3849
%C=[0.3828];
%C=[0.38266];
%C=[0.382709];
C=[0.382709087918741];
%Beta=0.55:0.01:0.65;
%Beta=[0.152];
%Beta=linspace(0.01,0.2,20);
%Beta=0.10:0.005:0.14; %Betamax 0.1667~0.1871
%Beta=[0.100051419253494];
%Beta=[0.100052559863590];
Beta=[0.100052559862974];
F=zeros(length(C),length(Beta));
X=[];
FX=[];
X1=[];
FX1=[];
X2=[];
FX2=[];
for i=1:length(C)
    for j=1:length(Beta)
        c=C(i);
        beta=Beta(j);
x = optimvar('x',1,9);
%Create the objective function as an expression of the optimization variable.

xy=[x(4),x(5),x(6),x(7),x(8),x(9)]'*[x(4),x(5),x(6),x(7),x(8),x(9)];
hxy=-xy.*log(xy)-(1-xy).*log(1-xy);
qb=1-x(3);
ehxy=[qb*x(1),qb*x(2),qb*(1-x(1)-x(2)),x(3)*x(1),x(3)*x(2),x(3)*(1-x(1)-x(2))]*hxy*[qb*x(1),qb*x(2),qb*(1-x(1)-x(2)),x(3)*x(1),x(3)*x(2),x(3)*(1-x(1)-x(2))]';

pi0=[x(4),x(5),x(6)]'*[x(4),x(5),x(6)]+[x(4)*(1-x(4)),x(5)*(1-x(5)),x(6)*(1-x(6))]'*[x(4)*(1-x(4)),x(5)*(1-x(5)),x(6)*(1-x(6))];
hpi0= -pi0.*log(pi0)-(1-pi0).*log(1-pi0);
pi1=[x(7),x(8),x(9)]'*[x(7),x(8),x(9)]+[x(7)*(1-x(7)),x(8)*(1-x(8)),x(9)*(1-x(9))]'*[x(7)*(1-x(7)),x(8)*(1-x(8)),x(9)*(1-x(9))];
hpi1= -pi1.*log(pi1)-(1-pi1).*log(1-pi1);
ehpi=[x(1),x(2),1-x(1)-x(2)]*(qb*hpi0+x(3)*hpi1)*[x(1),x(2),1-x(1)-x(2)]';

xx=[x(4),x(5),x(6),x(7),x(8),x(9)];
hx=-xx.*log(xx)-(1-xx).*log(1-xx);
ehx=[qb*x(1),qb*x(2),qb*(1-x(1)-x(2)),x(3)*x(1),x(3)*x(2),x(3)*(1-x(1)-x(2))]*hx';

obj = ((1-beta)*ehxy+beta*ehpi)/ehx;

%Create an optimization problem named prob having obj as the objective function.
prob = optimproblem('Objective',obj);

% options = optimoptions("fmincon",...
%     "Algorithm","interior-point",...
%     "EnableFeasibilityMode",true,...
%     'TolX', 1e-8, 'TolFun', 1e-8,...
%     'MaxIter',1000,...
%     "SubproblemAlgorithm","cg");

%options = optimoptions("fmincon",'Algorithm', 'interior-point',...
%    'TolX', 1e-8, 'TolFun', 1e-8);

options = optimoptions('fmincon','Algorithm','sqp', 'TolCon', 1e-12, 'TolX', 1e-12, 'OptimalityTolerance',1e-12,'MaxIterations', 4000,'FunctionTolerance',1e-8);

%Create the nonlinear constraint 
meancons =qb*(x(1)*x(4)+x(2)*x(5)+(1-x(1)-x(2))*x(6)) + x(3)*(x(1)*x(7)+x(2)*x(8)+(1-x(1)-x(2))*x(9)) >=1-c;
lb0 = x(4) >=0;
lb1 = x(5) >=0;
lb2 = x(6) >=0;
lb3 = x(7) >=0;
lb4 = x(8) >=0;
lb5 = x(9) >=0;
ub0 = x(4) <=1;
ub1 = x(5) <=1;
ub2 = x(6) <=1;
ub3 = x(7) <=1;
ub4 = x(8) <=1;
ub5 = x(9) <=1;
la1 = x(1) >=0;
ua1 = x(1) <=1;
la2 = x(2) >=0;
ua2 = x(2) <=1;
ca12 = x(1)+x(2) <=1;
lq = x(3) >=0;
uq = x(3) <=1;


%Include the nonlinear constraint in the problem. 
prob.Constraints.meancons =meancons;
prob.Constraints.lb0 = lb0;
prob.Constraints.lb1 = lb1;
prob.Constraints.lb2 = lb2;
prob.Constraints.lb3 = lb3;
prob.Constraints.lb4 = lb4;
prob.Constraints.lb5 = lb5;
prob.Constraints.ub0 = ub0;
prob.Constraints.ub1 = ub1;
prob.Constraints.ub2 = ub2;
prob.Constraints.ub3 = ub3;
prob.Constraints.ub4 = ub4;
prob.Constraints.ub5 = ub5;
prob.Constraints.la1 = la1;
prob.Constraints.ua1 = ua1;
prob.Constraints.la2 = la2;
prob.Constraints.ua2 = ua2;
prob.Constraints.ca12 = ca12;
prob.Constraints.lq = lq;
prob.Constraints.uq = uq;

%Review the problem.
%show(prob)
        F(i,j)=2;
        for it=1:4000% 100iter gives spurious with prob 1/15
        x0.x = rand(1,9);
        a1=x0.x(1);
        a2=x0.x(2);
        a3=rand;
        x0.x(1)=a1/(a1+a2+a3);
        x0.x(2)=a2/(a1+a2+a3);
        %x0.x=local22(:,1)';
            try
                [sol,fval,exitflag,output] = solve(prob,x0,'Options', options);
            catch exception
                fval=3;
            end
        F(i,j)=min(fval,F(i,j));  
            if fval<1
                X=[X;sol.x];
                FX=[FX;fval];    
            end
            %if fval<1.00008059
            if fval<1.000002
                X1=[X1;sol.x];
                FX1=[FX1;fval];    
            end
            if fval<1.002
                X2=[X2;sol.x];
                FX2=[FX2;fval];    
            end
        end
    end
end
%[CC,BBeta] = meshgrid(C,Beta);
%mesh(CC',BBeta',F)
%plot(Beta,F(5,:))
%plot(C,max(F'))
fprintf('%.6f\n', F');


%Beta=linspace(0.01,0.2,20);
%F=zeros(length(C),length(Beta));
% 0.9994    0.9998    1.0001    1.0004    1.0007    1.0010    1.0013    1.0016    1.0019    1.0022    1.0025
%     0.9990    0.9993    0.9996    0.9999    1.0002    1.0005    1.0008    1.0011    1.0014    1.0017    1.0021
% 
%   Columns 12 through 20
% 
%     1.0028    1.0031    1.0035    1.0038    1.0041    1.0044    1.0047    1.0039    0.9997
%     1.0024    1.0027    1.0030    1.0033    1.0036    1.0039    0.9995    1.0045    1.0048

% C=0.3849
% [Beta',F']=
% 0.1000    0.9984
%     0.1500    1.0000
%     0.2000    1.0016
%     0.2500    0.9987
%     0.3000    1.0048
%     0.3500    0.9937
%     0.4000    0.9928
%     0.4500    0.9918
%     0.5000         0
%another trial with 200 it:
% 0.1000    0.9984
%     0.1500    0.9963
%     0.2000    0.9959
%     0.2500    0.9952
%     0.3000    0.9945
%     0.3500    0.9937
%     0.4000    0.7462
%     0.4500    0.9918
%     0.5000    0.9908

%C=[0.3828]; Beta=0.1:0.02:0.18; it=200 start 9:40 finish before 10:17
%F = 1.0017    1.0024    1.0030    1.0036    0.9995
% anothor trial start 10:19 end before 10:45
%1.0017    1.0024    0.9997    1.0036    XXX1.0042
%Beta=0.02:0.02:0.14 start10:50 end 11:08  2800-18min
%0.9993    0.9999    1.0005    1.0011    1.0017    1.0024    0.9997
%Beta=0.03:0.02:0.13; start 11:12 end 11:27
%0.9996    1.0002    1.0008    1.0014    1.0021    1.0027
%Beta=0.10:0.005:0.14 start 11:32 stopped manually
%F= 1.001747417995477   1.001902251814894   1.002057085603230
%then I change TolX and TolFun to 1e-8 (before they were 1e-6 and 1e-4)
%%%%%%  repeat 11:45
% 1.001747
% 1.001902
% 1.002057
% 0.999837
% 1.002367XX
% 1.002522XX
% 1.002676XX
% 1.002831XX
% 1.002986XX
%change it=500, Beta=0.10:0.005:0.12
% 1.001747
% 1.001902
% 0.999846
% 1.002212
% 1.002367
% 1.002522
% 1.002676
% 0.999767
% 0.999742
%repeat 1:32 to check why solutions <1
% 1.001747
% 1.001902
% 1.002057
% 1.002212
% 1.002367
% 1.002522
% 1.002676
% 0.999767
% 1.002986

% local min: the first entry in variable X
% ehxy/ehx
% ans =
%    1.000403935820589
% ehpi/ehx
% ans =
%    0.995683828685120

%beta, c=optimal, it=4100; found opt fval=1.000080506484608
%another 5000, iter 3 got to this minimum
% FX1 =
%    1.000080506487452
%    1.000080506484608
%    1.000080506485429
% X1'=0.446836396219786   0.106323249682956   0.446836986215441
%    0.446838641430491   0.313727993778602   0.446838047561532
%    0.999777263739876   0.999777359083513   0.999777264067821
%    0.688561470492352   0.681161746436530   0.688561398977144
%    0.688561530153092   0.687564249941912   0.688561395745462
%    0.681176028472391   0.689107372190856   0.681176845590311
%    0.690770614944445   0.000000474889695   0.690770619059044
%    0.690770613978640   0.690769155304157   0.690770615906549
%    0.000000474880894   0.690769374604828   0.000000474879054

%new test C=[0.382709];1000iter
%27000 iter F=1.000001127354140 good
% spurious: fval->1.001892196451269
%    0.309095796575199 
%    0.233160365963954
%    0.441702614567072
%    0.617290151054174
%    0.617290460375643
%    0.617291684898728
%    0.617290485767789
%    0.617291739590452
%    0.617291428387402
% another 27000 iter
% FX1 =
%    1.000001127380566
%    1.000001127379323
%    1.000001127385291
%    1.000001127381497
%    1.000001127392111
%    1.000001127380090
%    1.000001127380953
%    1.000001127379062
%    1.000001127380904
%X1(end,;)
%  0.106394141274807
%    0.313672719256056
%    0.999777675096006
%    0.681163865560585
%    0.687550977216031
%    0.689100757936282
%    0.000000474546332
%    0.690769049714938
%    0.690769402295826
% change tol to e-16, 3000 it:
% FX1 =
%    1.000000984954297
%    1.000001002724186
% X1'=
%    0.446806368953906   0.106394644449426
%    0.446796926885473   0.322461681826042
%    0.999777588074144   0.000219839965392
%    0.688557319723004   0.000000475366016
%    0.688557311325121   0.690767541102705
%    0.681169552481399   0.690771447563901
%    0.690771157838579   0.677736902398077
%    0.690771150438346   0.686825745520445
%    0.000000474531305   0.688318330462258

%run 5000 iterations
% X1'=
%    0.446801401300231   0.106395079138204
%    0.446801910804794   0.313678345457482
%    0.999777587867945   0.000222324407373
%    0.688557177195678   0.000000474540767
%    0.688560182721519   0.690769672017182
%    0.681166548407982   0.690770033702867
%    0.690771141600706   0.681155280634292
%    0.690771141048416   0.687546181946635
%    0.000000474531335   0.689093894825876
% FX1= 1.000000984954580   1.000000984974595
%Tol changed to 1e-9,
%    0.106394971482438
%    0.313667201426177
%    0.000222318425183
%    0.000000474541749
%    0.690769600268406
%    0.690769942341982
%    0.681159342081384
%    0.687556807037114
%    0.689099322652099

%switch to sqp, 4000it
%sqp-legacy 400 it

