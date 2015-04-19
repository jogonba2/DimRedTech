function [m,W] = pca(X,k)
	# Calcular la media de los vectores #
    m                = (1/size(X,2)) * sum(X,2);
    # Restar a cada vector la media #
    t                = X .- repmat(m,1,size(X,2));
    # Calcular la matriz de covarianza #
    matrizCovarianza = (1/size(t,2)) * (t*t');
    # Calcular los eigenvalores y eigenvectores de la matriz de covarianza #
    [eigenVectores,eigenValores] = eig(matrizCovarianza);
    # Ordenar eigenvectores por mayor eigenvalor #
    [ordered,perm]   = sort(diag(eigenValores),'descend');
    eigenVectores    = eigenVectores(:,perm);
    # Definir la matriz W como los k primeros eigenvectores #
    W                = eigenVectores(:,[1:k]);
end

# Para imprimir el eigenvector K -> 
#>> Cargar el fichero en X, luego X=X'
#>> [media,eigenvectores] = pca(X,k)
#>> x=eigenvectores(:,k);
#>> xr=reshape(x,16,16);
#>> imshow(xr,[min(min(xr)),max(max(xr))]) // Imprimir

# Ordenar eigenvectores por eigenvalor -> https://ddcampayo.wordpress.com/2015/02/20/sorting-eivenvectors-in-octave/ #
