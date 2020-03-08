clear all; close all; clc;

ff = 0.0567;

nn = 10;

weight_vec = weight_vec(ff,nn);

sum(weight_vec)

figure(1); plot(weight_vec);