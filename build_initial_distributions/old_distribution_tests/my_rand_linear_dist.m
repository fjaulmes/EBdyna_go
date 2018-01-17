%function to return a rand with a linear decreasing proba dist (finf>fsup)
function invF=my_rand_linear_dist(finf,fsup,xinf,xsup)

if fsup~=0
    beta=finf/fsup;
else
    fsup=0.1*finf;
    beta=finf/fsup;
end
if beta~=-1
    a=(1-beta)/(1+beta);
else
    a=0;
end
    
b=1-a;
c=0;

if a~=0
    cr=-my_rand(1);
    delta=b^2-4*a*cr;
    
    invF=(-b+sqrt(delta))/(2*a);
    invF=invF*(xsup-xinf)+xinf;
else
    cr=my_rand(1);
    invF=cr*(xsup-xinf)+xinf;
end
