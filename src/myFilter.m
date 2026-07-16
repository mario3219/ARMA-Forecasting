% na = AR order, nc = MA order

function arma_model = myFilter(data, orders, noLags, titleStr)

arma_model = armax(data, [orders(1) orders(2)]);
e_hat = filter(arma_model.A, arma_model.C, data);
plotACFnPACF(e_hat, noLags, titleStr);