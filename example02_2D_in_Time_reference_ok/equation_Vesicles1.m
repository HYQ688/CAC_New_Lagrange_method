%% Epitaxy
clear 
clc
syms x y epsilon t Lx Ly alpha beta gamma Re lambda M eta2

phi = exp(-t)*sin(2*x)*cos(2*y);
u =  sin(2*y).*sin(x).^2.*sin(t);
v = -sin(2*x).*sin(y).^2.*sin(t); 
p  = cos(x).*sin(y).*sin(t);

div_u = simplify(diff(u,x,1)+ diff(v,y,1));
lap_phi = diff(phi,x,2) + diff(phi,y,2);
lap2_phi = diff(lap_phi,x,2) + diff(lap_phi,y,2);
lap_u = diff(u,x,2) + diff(u,y,2);
lap_v = diff(v,x,2) + diff(v,y,2);
grad_square = diff(phi,x,1).^2 + diff(phi,y,1).^2;

phi1 = diff(phi,x,1); phi2 = diff(phi,y,1);
u1 = diff(u,x,1); u2 = diff(u,y,1);
v1 = diff(v,x,1); v2 = diff(v,y,1);
p1 = diff(p,x,1); p2 = diff(p,y,1);

advection   = diff(u.*phi,x,1) + diff(v.*phi,y,1);
convective1 = u.*u1 + v.*u2;
convective2 = u.*v1 + v.*v2;

F     = (phi.^2-1).^2/4;
f     = phi.^3 - phi;
f_der = 3*phi.^2 - 1;

omega = -lap_phi +1./epsilon.^2.*f;
lap_omega = diff(omega,x,2) + diff(omega,y,2);

% A_phi = int(int(phi+1,x,0,2*pi),y,0,2*pi);
% B_phi = int(int(epsilon./2.*grad_square+1./epsilon.*F,x,0,2*pi),y,0,2*pi);

mu = lambda*epsilon*( lap2_phi+6./epsilon*( grad_square*phi - diff(phi^2*phi1,x,1) - diff(phi^2*phi2,y,1) ) + 1./epsilon^4*f_der*f + 2./epsilon^2*lap_phi ) -lambda*eta2*(-epsilon*lap_phi+1./epsilon*f);

tension1    =  mu.*diff(phi,x,1);
tension2    =  mu.*diff(phi,y,1);

mubar = mu - 1./(2*pi*2*pi)*int(int(mu,x,0,2*pi),y,0,2*pi);

phit = simplify( diff(phi,t));
ut = simplify( diff(u,t));
vt = simplify( diff(v,t));

f1 = simplify( phit + advection + gamma*mubar);
f2 = simplify( ut  + convective1 + p1 - lap_u/Re - tension1);
f3 = simplify( vt  + convective2 + p2 - lap_v/Re - tension2);


f1 = char(f1); f1 = strrep(f1,'*','.*'); f1 = strrep(f1,'/','./');
f1 = strrep(f1,'^','.^')

f2 = char(f2); f2 = strrep(f2,'*','.*'); f2 = strrep(f2,'/','./');
f2 = strrep(f2,'^','.^')

f3 = char(f3); f3 = strrep(f3,'*','.*'); f3 = strrep(f3,'/','./');
f3 = strrep(f3,'^','.^')

phi 

phit = char(phit);
phit = strrep(phit,'*','.*');
phit = strrep(phit,'/','./');
phit = strrep(phit,'^','.^')

u
ut = char(ut);
ut = strrep(ut,'*','.*');
ut = strrep(ut,'/','./');
ut = strrep(ut,'^','.^')

v
vt = char(vt);
vt = strrep(vt,'*','.*');
vt = strrep(vt,'/','./');
vt = strrep(vt,'^','.^')

p
p1
p2
