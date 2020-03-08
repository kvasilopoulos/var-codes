function weight_vec=weight_vec(ff,nn)

d_vec   =   ones(1,2*nn);

for kk=1:2*nn
    d_vec(1,kk)=d(ff,kk/2);
end;

weight_H_vec = NaN(1,2*nn);

for kk=1:2*nn
    weight_H_vec(1,kk)=(d_vec(1,2*nn)+0.5*ff*sum(d_vec(1,kk:2*nn)))/(0.5*sum(d_vec(1,1:2*nn)));
end;

weight_vec = NaN(1,12*nn);

for kk=1:12*nn
    HH = floor((kk-1)/6)+1;
    weight_vec(1,kk)=weight_H_vec(1,HH)/12;
end;


