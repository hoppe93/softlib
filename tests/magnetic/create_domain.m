clc; clear;

Rmajor = 3;
Rminor = 1;

figure('Position', [500, 400, 600, 1400]);
lim =[Rmajor-Rminor Rmajor+Rminor -2*Rminor 2*Rminor];
axis(lim);

waiting = true;
X = [];
Y = [];

while waiting
    [x,y,b] = ginput(1);
    
    if b == 1
        X = [X x];
        Y = [Y y];
    elseif b == 2
        waiting = false;
    elseif numel(X) > 0
        X = X(1:end-1);
        Y = Y(1:end-1);
    end
    
    clf;
    hold on
    plot(X, Y, 'linewidth', 2);
    plot(X, Y, 'r.');
    axis(lim);
end

% Return to final point
X = [X X(1)];
Y = [Y Y(1)];
    
clf;
hold on
plot(X, Y, 'linewidth', 2);
plot(X, Y, 'r.');
axis(lim);

% Select some points inside the domain
disp('Select points inside the domain');
iX = [];
iY = [];
waiting = true;
while waiting
    [x,y,b] = ginput(1);
    
    if b == 1
        iX = [iX x];
        iY = [iY y];
    elseif b == 2
        waiting = false;
    elseif numel(iX) > 0
        iX = iX(1:end-1);
        iY = iY(1:end-1);
    end
    
    clf;
    hold on
    plot(X, Y, 'linewidth', 2);
    plot(X, Y, 'r.');
    plot(iX, iY, 'rx', 'linewidth', 2);
    plot(iX, iY, 'g--', 'linewidth', 2);
    axis(lim);
end

% Select points outside domain
disp('Select points outside the domain');
oX = [];
oY = [];
waiting = true;
while waiting
    [x,y,b] = ginput(1);
    
    if b == 1
        oX = [oX x];
        oY = [oY y];
    elseif b == 2
        waiting = false;
    elseif numel(oX) > 0
        oX = oX(1:end-1);
        oY = oY(1:end-1);
    end
    
    clf;
    hold on
    plot(X, Y, 'linewidth', 2);
    plot(X, Y, 'r.');
    plot(oX, oY, 'ro', 'linewidth', 2);
    plot(oX, oY, 'g--', 'linewidth', 2);
    axis(lim);
end

% Generate C arrays
s = ['const unsigned int NWALL=',num2str(numel(X)),';',10,'slibreal_t test_domain[2][NWALL] = {',10];
is= ['const unsigned int NINPOINTS=',num2str(numel(iX)),';',10,'slibreal_t test_domain_inside[2][NINPOINTS] = {',10];
os= ['const unsigned int NOUTPOINTS=',num2str(numel(oX)),';',10,'slibreal_t test_domain_outside[2][NOUTPOINTS] = {',10];
x = '    {';
y = '    {';
ix= '    {';
iy= '    {';
ox= '    {';
oy= '    {';
for i=1:numel(X)
    x = [x, num2str(X(i)),','];
    y = [y, num2str(Y(i)),','];
end
for i=1:numel(iX)
    ix = [ix, num2str(iX(i)),','];
    iy = [iy, num2str(iY(i)),','];
end
for i=1:numel(oX)
    ox = [ox, num2str(oX(i)),','];
    oy = [oy, num2str(oY(i)),','];
end
s = [s, x(1:end-1), '},', 10, y(1:end-1), '}', 10, '};'];
is= [is,ix(1:end-1),'},', 10,iy(1:end-1), '}', 10, '};'];
os= [os,ox(1:end-1),'},', 10,oy(1:end-1), '}', 10, '};'];

disp(s);
disp(is);
disp(os);
