function weight=weight(ff,nn)

d_vec   =   ones(1,2*nn);

for kk=1:2*nn
    d_vec(1,kk)=d(ff,kk/2);
end;

weight=1/(1-d_vec(1,2*nn)+0.5*sum(d_vec(1,1:2*nn)));



