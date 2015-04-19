function [W] = lda(X,xl,k)
# Media de vectores #
# xl en filas y X en columnas #
m  = (1/columns(X)) * sum(X,2);
Sw = zeros(rows(X));
Sb = zeros(rows(X));
for i=1:10
    indices = find(xl==i-1);
    # Calcular media de la clase i #
    mc      = (1/rows(indices)) * sum(X(:,indices),2);
    # Calcular matriz covarianza de la clase i #
    t       = X(:,indices) .- repmat(mc,1,columns(X(:,indices)));
    matrizCovarianza = (1/columns(t)) * (t*t');
    # Calcular matriz intra clases Sw #
    Sw      = Sw + matrizCovarianza;
    # Calcular matriz entre clases Sb #
    aux     = mc .- m;
    Sb      = Sb + rows(indices) * (aux*aux');
endfor
# Eigenvectores generalizados de Sb y Sw #
[eigenVectores,eigenValores] = eig(Sb,Sw);
# Ordenar eigenvectores por mayor eigenvalor #
[ordered,perm]   = sort(diag(eigenValores),'descend');
eigenVectores    = eigenVectores(:,perm);
# Definir la matriz W como los k primeros eigenvectores #
W                = eigenVectores(:,[1:k]);
end

# Para mostrarlo #
# load trlabels.mat
# load trdata.mat
# X = X';
# W = lda(X,xl,k);
# x = W(:,y);
# xr = reshape(x,16,16);
# imshow(xr,[min(min(xr)),max(max(xr))])
