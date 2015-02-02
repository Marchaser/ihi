function [u,du_dc,du_dl,du_dh] = U(c,l,h,lambda,eta,rho,nu)
cl_comp = c.^eta .* l.^eta;
dcl_comp_dc = eta .* cl_comp ./ c;
dcl_comp_dl = (1-eta) .* cl_comp ./ l;
u_base = lambda*cl_comp.^rho + (1-lambda)*h.^rho;
u_base_power = u_base.^((1-nu)/rho);
u_base_power_minus_1 = u_base_power ./ u_base;
u = u_base_power / (1-nu);
du_dh = (1-lambda) * h.^(rho-1) .* u_base_power_minus_1;
du_dcl_comp = lambda * u_base_power_minus_1 .* cl_comp.^(rho-1);
du_dc = du_dcl_comp .* dcl_comp_dc;
du_dl = du_dcl_comp .* dcl_comp_dl;
end