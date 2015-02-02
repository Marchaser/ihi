function [u,du_dc,du_dl,du_dh] = UHJ(c,l,h,lambda,zeta,rho,nu)
cl_comp = c.^zeta .* l.^zeta;
dcl_comp_dc = zeta .* cl_comp ./ c;
dcl_comp_dl = (1-zeta) .* cl_comp ./ l;
u = cl_comp.^(1-nu)/(1-nu) + lambda*h.^(1-rho)/(1-rho);
du_dh = lambda*h.^(-rho);
du_dcl_comp = cl_comp.^(-nu);
du_dc = du_dcl_comp .* dcl_comp_dc;
du_dl = du_dcl_comp .* dcl_comp_dl;
end