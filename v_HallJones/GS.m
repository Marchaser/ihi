function tax = GS(y, tau0, tau1, tau2)
bracket1 = y.^(-tau1) + tau2;
bracket2 = y - bracket1.^(-1/tau1);
tax = tau0 * bracket2;
end