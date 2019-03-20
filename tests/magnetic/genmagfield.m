clc; clear;

npoints = 25;
Rm = 0.68;
B0 = 5;
rminor = 0.22;

% q-profile parameters (a1 : c, a2 : l, a3 : q, a4 : e)
a1 = 0.5;
a2 = 2;
a3 = 2;
a4 = log(2);

% q-profiles
qc = @(r) a1;
ql = @(r) a2*r + 1;
qq = @(r) a3*r.^2 + 1;
qe = @(r) exp(a4*r);

% r-derivatives of q-profiles (times r)
qpc= @(r) 0;
qpl= @(r) a2*r;
qpq= @(r) 2*a3*r.^2;
qpe= @(r) r.*exp(r*a4) * a4;

%% Generate .cpp-files

% Buffers
sc = [];
sl = [];
sq = [];
se = [];
% Precision in number-to-string conversion
ndigits = 15;

ns = @(s) num2str(s, ndigits);

% Generate random points in 3D space
for i=1:npoints
    minr = rminor * rand;
    pola = 2*pi * rand;
    tora = 2*pi * rand;
    
    x = [(Rm-minr*cos(pola))*cos(tora), (Rm-minr*cos(pola))*sin(tora), minr*sin(pola)];
    
    [Bc, cBabs, gradBc, curlBc] = magnetic_field(x, B0, Rm, -1, +1, qc(minr/rminor), qpc(minr/rminor));
    [Bl, lBabs, gradBl, curlBl] = magnetic_field(x, B0, Rm, -1, +1, ql(minr/rminor), qpl(minr/rminor));
    [Bq, qBabs, gradBq, curlBq] = magnetic_field(x, B0, Rm, -1, +1, qq(minr/rminor), qpq(minr/rminor));
    [Be, eBabs, gradBe, curlBe] = magnetic_field(x, B0, Rm, -1, +1, qe(minr/rminor), qpe(minr/rminor));
    
    c = [x, Bc', gradBc', curlBc'];
    l = [x, Bl', gradBl', curlBl'];
    q = [x, Bq', gradBq', curlBq'];
    e = [x, Be', gradBe', curlBe'];
    
    sc = [sc, '\t{'];
    for j=1:numel(c)
        sc = [sc,ns(c(j)),','];
    end
    sc(end) = '}';
    
    sl = [sl, '\t{'];
    for j=1:numel(l)
        sl = [sl,ns(l(j)),','];
    end
    sl(end) = '}';
    
    sq = [sq, '\t{'];
    for j=1:numel(q)
        sq = [sq,ns(q(j)),','];
    end
    sq(end) = '}';
    
    se = [se, '\t{'];
    for j=1:numel(c)
        se = [se,ns(e(j)),','];
    end
    se(end) = '}';
    
    if i < npoints
        sc = [sc,',\n'];
        sl = [sl,',\n'];
        sq = [sq,',\n'];
        se = [se,',\n'];
    else
        sc = [sc,'\n'];
        sl = [sl,'\n'];
        sq = [sq,'\n'];
        se = [se,'\n'];
    end
end

% Generate header file
fid = fopen('magfield_points.h', 'w');

fwrite(fid, ['/* This file was auto-generated */',10]);
fwrite(fid, ['#ifndef _MAGFIELD_POINTS_H',10,'#define _MAGFIELD_POINTS_H',10,10]);
fwrite(fid, ['#include <softlib/config.h>',10,10]);
fwrite(fid, ['#define MAGNETIC_FIELD_TEST_GUARANTEED_PRECISION (max(100*',num2str(10^(-ndigits)),',20*REAL_EPSILON))',10,10]);
fwrite(fid, ['const unsigned int MAGNETIC_FIELD_TEST_NPOINTS=',num2str(npoints),';',10]);
fwrite(fid, ['const slibreal_t magnetic_field_test_data_B0      =',ns(B0),',',10]);
fwrite(fid, ['                 magnetic_field_test_data_Rm      =',ns(Rm),',',10]);
fwrite(fid, ['                 magnetic_field_test_data_rminor  =',ns(rminor),',',10]);
fwrite(fid, ['                 magnetic_field_test_data_qa_const=',ns(a1),',',10]);
fwrite(fid, ['                 magnetic_field_test_data_qa_lin  =',ns(a2),',',10]);
fwrite(fid, ['                 magnetic_field_test_data_qa_quad =',ns(a3),',',10]);
fwrite(fid, ['                 magnetic_field_test_data_qa_exp  =',ns(a4),';',10]);
fwrite(fid, ['extern const slibreal_t magnetic_field_test_data_const[',num2str(npoints),'][12];',10]);
fwrite(fid, ['extern const slibreal_t magnetic_field_test_data_linear[',num2str(npoints),'][12];',10]);
fwrite(fid, ['extern const slibreal_t magnetic_field_test_data_quadratic[',num2str(npoints),'][12];',10]);
fwrite(fid, ['extern const slibreal_t magnetic_field_test_data_exponential[',num2str(npoints),'][12];',10,10]);

fwrite(fid, '#endif/*_MAGFIELD_POINTS_H*/');
fclose(fid);

% Generate table
fid = fopen('magfield_points.cpp', 'w');

fwrite(fid, ['/* This file was auto-generated */',10,10]);
fwrite(fid, ['#include <softlib/config.h>',10,'#include "magfield_points.h"',10,10]);
%fwrite(fid, ['const unsigned int MAGNETIC_FIELD_TEST_NPOINTS=',num2str(npoints),';',10,10]);

fwrite(fid, ['/* Format = X, B(X), grad|B(X)|, curl B(X) */',10,10]);

fwrite(fid, ['const slibreal_t magnetic_field_test_data_const[',num2str(npoints),'][12] = {',10]);
fprintf(fid, sc);
fwrite(fid, ['};',10,10]);

fwrite(fid, ['const slibreal_t magnetic_field_test_data_linear[',num2str(npoints),'][12] = {',10]);
fprintf(fid, sl);
fwrite(fid, ['};',10,10]);

fwrite(fid, ['const slibreal_t magnetic_field_test_data_quadratic[',num2str(npoints),'][12] = {',10]);
fprintf(fid, sq);
fwrite(fid, ['};',10,10]);

fwrite(fid, ['const slibreal_t magnetic_field_test_data_exponential[',num2str(npoints),'][12] = {',10]);
fprintf(fid, se);
fwrite(fid, ['};',10,10]);

fclose(fid);

%% Generate a SOFT magnetic equilibrium file

nr = 100;
nz = 100;
r = linspace(-rminor, rminor, nr) + Rm;
z = linspace(-rminor, rminor, nz);
name = 'Circular magnetic field - for benchmarking SOFT';
desc = ['Circular magnetic field with B0 = ',num2str(B0,2),', Rm = ',num2str(Rm,2),', q = ',num2str(a1,2)];

% Generate wall
t = linspace(0, 2*pi);
wall = [Rm+rminor*cos(t); rminor*sin(t)];

Br = zeros(nr, nz);
Bphi = zeros(nr, nz);
Bz = zeros(nr, nz);
maxis = [Rm 0];

for i=1:nr
    for j=1:nz
        minr = r(i) - Rm;
        B = magnetic_field([r(i), 0, z(j)], B0, Rm, -1, +1, qc(minr/rminor), qpc(minr/rminor));
        
        Br(i,j)   = B(1);
        Bphi(i,j) = B(2);
        Bz(i,j)   = B(3);
    end
end

save('circular-benchmark.mat', '-v7.3', 'Bphi', 'Br', 'Bz', 'desc', 'maxis', 'name', 'r', 'wall', 'z');
