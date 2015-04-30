#!/usr/bin/octave -qf
if(nargin!=9)
    printf("Usage: pcalda.m <type> <trdata> <trlabels> <tsdata> <tslabels> <mink> <stepk> <maxk> <result_file> \n");
    exit(1);
end

# Argumentos #
arg_list=argv();
type     = arg_list{1};
trdata   = arg_list{2};
trlabs   = arg_list{3};
tsdata   = arg_list{4};
tslabs   = arg_list{5};
mink     = str2num(arg_list{6});
stepk    = str2num(arg_list{7});
maxk     = str2num(arg_list{8});
nameFile = arg_list{9};
# Carga de archivos #
load(trdata);
load(trlabs);
load(tsdata);
load(tslabs);
X = X';
Y = Y';

# Gestión de script #

## Espacio Original (Solo se calculará una vez) ##
if(type=="ORI")
	printf("Error en el espacio original : %f\n",knn(X,xl,Y,yl,1));
end

fd = fopen(nameFile,"w");
## PCA ##
if(type=="PCA")
	fwrite(fd,"# This file has PCA results #\r\n");
	fwrite(fd,"# Dimensions\t%ErrorPCA\r\n");
    while mink<=maxk
		[avg,W] = pca(X,mink);
		proyX   = W'*X;
		proyY   = W'*Y;
		fwrite(fd,num2str(mink));
		fwrite(fd,"\t\t");
		fwrite(fd,num2str(knn(proyX,xl,proyY,yl,1)));
		fwrite(fd,"\r\n");
		mink = mink+stepk;
	endwhile
	fclose(fd);
end

## LDA (Se asumen parámetros correctos - k <= C-1 -) ##
if(type=="LDA")
	fwrite(fd,"# This file has LDA results #\r\n");
	fwrite(fd,"# Dimensions\t%ErrorLDA\r\n");
	while mink<=maxk
		W		= lda(X,xl,mink);
		proyX   = W'*X;
		proyY	= W'*Y;
		fwrite(fd,num2str(mink));
		fwrite(fd,"\t\t");
		fwrite(fd,num2str(knn(proyX,xl,proyY,yl,1)));
		fwrite(fd,"\r\n");
		mink = mink+stepk;
	endwhile
	fclose(fd);
end

## PCA y LDA separados generando la gráfica ##
if(type=="PLA")
	fwrite(fd,"# This file has PCA and LDA results separated#\r\n");
	fwrite(fd,"# Dimensions\t%ErrorPCA\t%ErrorLDA\r\n");
	while mink<=maxk
		[m,WPca]  = pca(X,mink);
		WLda	  = lda(X,xl,mink);
		proyPcaX  = WPca' * X;
		proyPcaY  = WPca' * Y;
		proyLdaX  = WLda' * X;
		proyLdaY  = WLda' * Y;
		fwrite(fd,num2str(mink));
		fwrite(fd,"\t\t");
		fwrite(fd,num2str(knn(proyPcaX,xl,proyPcaY,yl,1)));
		fwrite(fd,"\t\t");
		fwrite(fd,num2str(knn(proyLdaX,xl,proyLdaY,yl,1)));
		fwrite(fd,"\r\n");
		mink = mink+stepk;
	endwhile
	fclose(fd);
	system("gnuplot PCALDA.gp");
end	



# PCA y LDA combinados método 2. Usar también en la competición #
if(type=="PLM")
    minPCA = 1;
    maxPCA = 256; # Ajustar dimensiones PCA #
    jump   = 10;
    minOneDimPCA    = 0;
    minTwoDimPCA    = 0;
    minOneErrorPCA  = 100;
    minTwoErrorPCA  = 100;
    kk = 1;
    # Calcular las 2 dimensiones de PCA que dan menor error (entre 0 y 256)
    while minPCA<maxPCA
		[avg,W] = pca(X,minPCA);
		proyX   = W'*X;
		proyY   = W'*Y;
		err     = knn(proyX,xl,proyY,yl,kk);
		if(err<minOneErrorPCA)
			# Si el error es menor que el primer valor del intervalo, el menor anterior pasa a ser 
			# el 2º menor y se actualiza el 1º menor
			minTwoDimPCA = minOneDimPCA;
		    minTwoErrorPCA = minOneErrorPCA;
		    minOneDimPCA = minPCA;
		    minOneErrorPCA = err;		    
		else
		    if(err<minTwoErrorPCA)
		        minTwoDimPCA = minPCA;
		        minTwoErrorPCA = err;
		    end
		end	 
		if(minPCA+jump<maxPCA) minPCA = minPCA + jump;
		else minPCA = minPCA + 1;
		end
		disp(minPCA);
    endwhile
    disp("Menor error 1: "),disp(minOneErrorPCA),disp(" en la dimension: "),disp(minOneDimPCA);
    disp("Menor error 2: "),disp(minTwoErrorPCA),disp(" en la dimension: "),disp(minTwoDimPCA);
    # Para cada valor entre las 2 dimensiones de menor error de PCA, calcular su error con proyección LDA #
    dimPCA = 1;
    dimLDA = 1;
    menorError = 100;
    while(minOneDimPCA<minTwoDimPCA)
        for j=1:9 # Probar con nº de clases dimensiones de LDA (Ajustar dimensiones LDA AQUI)#
	        [avg,WPCA]  = pca(X,minOneDimPCA);
	        proyXPCA = WPCA'*X;
	        proyYPCA = WPCA'*Y;
	        WLDA     = lda(proyXPCA,xl,j);
	        proyXLDA = WLDA'*proyXPCA;
	        proyYLDA = WLDA'*proyYPCA;
	        err      = knn(proyXLDA,xl,proyYLDA,yl,kk);
	        #disp("Calculado"),disp(minOneDimPCA),disp(j),disp("con error: "),disp(err);
	        if(err<menorError)
	            dimPCA = minOneDimPCA;
	            dimLDA = j;
	            menorError = err;
	            #disp("Menor error con: "),disp(dimPCA),disp(dimLDA),disp(" : "),disp(err)
	        end   
	    endfor
	    minOneDimPCA = minOneDimPCA + 1;
	endwhile
    disp("Menor error con: "),disp(dimPCA),disp(dimLDA),disp(" : "),disp(menorError)
end

## PCA Y LDA combinados (No son necesarios los 4 ultimos parametros), configurar a mano en la competicion ##
## Primero se aplica PCA proyectando hasta nº dimensiones de la imagen, después se aplica LDA sobre las 5 dimensiones que con PCA 
## den menor error ##
if(type=="PLC")
	minPCA     = 1;   # Ajustar dimension minima de PCA AQUI #
	maxPCA     = 256; # Ajustar dimension maxima de PCA AQUI #
	lenSet     = 10;
	menoresErrorPCA = [];
	menoresDimPCA   = [];
	kk         = 1; # Configurar numero de vecinos del clasificador #
	while minPCA<maxPCA
		[avg,W] = pca(X,minPCA);
		proyX   = W'*X;
		proyY   = W'*Y;
		err     = knn(proyX,xl,proyY,yl,kk);
		if(minPCA<=lenSet) # Coger las lenSet dimensiones que den menor error con PCA #
		    menoresErrorPCA(end+1) = err;
		    menoresDimPCA(end+1) = minPCA;
		else
		    if(err<max(menoresErrorPCA))
		        i = find(menoresErrorPCA==max(menoresErrorPCA))(1);
		        menoresErrorPCA(i) = err;
		        menoresDimPCA(i) = minPCA;
		    end
		end
		minPCA = minPCA + 1;
	#	disp("Con error: "),disp(menoresErrorPCA),disp(" dimension: "),disp(menoresDimPCA)
	endwhile
	disp("Menores con PCA ya calculados: "),disp(menoresDimPCA),disp(menoresErrorPCA)
	# Calcular con los 5 menores de PCA el error con LDA y coger el menor #
	x = 1;
	j = 1;
	#menoresDimPCA = [10,20,30,40,50];
	menorError = 100;
	dimPCA = 0;
	dimLDA = 0;
	for x=1:lenSet # Con las lenSet dimensiones de PCA que dan menor error #
	    for j=1:9 # Probar con nº de clases dimensiones de LDA (Ajustar dimensiones LDA AQUI)#
	        [avg,WPCA]  = pca(X,menoresDimPCA(x));
	        proyXPCA = WPCA'*X;
	        proyYPCA = WPCA'*Y;
	        WLDA     = lda(proyXPCA,xl,j);
	        proyXLDA = WLDA'*proyXPCA;
	        proyYLDA = WLDA'*proyYPCA;
	        err      = knn(proyXLDA,xl,proyYLDA,yl,kk);
	        #disp("Calculado"),disp(x),disp(j),disp("con error: "),disp(err);
	        if(err<menorError)
	            dimPCA = x;
	            dimLDA = j;
	            menorError = err;
	            #disp("Menor error con: "),disp(dimPCA),disp(dimLDA),disp(" : "),disp(err)
	        end   
	     endfor
	endfor
	disp("Menor error con: "),disp(menoresDimPCA(dimPCA)),disp(dimLDA),disp(" : "),disp(menorError)
	fclose(fd); 
end	
	
## PCA y Original generando la gráfica ##
if(type=="PCO")
	fwrite(fd,"# This file has PCA and Original results #\r\n");
	fwrite(fd,"# Dimensions\t%ErrorPCA\t%ErrorOriginal\r\n");
	while mink<=maxk
		[m,WPca]  = pca(X,mink);
		proyPcaX  = WPca' * X;
		proyPcaY  = WPca' * Y;
		fwrite(fd,num2str(mink));
		fwrite(fd,"\t\t");
		fwrite(fd,num2str(knn(proyPcaX,xl,proyPcaY,yl,1)));
		fwrite(fd,"\t\t");
		fwrite(fd,num2str(knn(X,xl,Y,yl,1)));
		fwrite(fd,"\r\n");
		mink = mink+stepk;
	endwhile
	fclose(fd);
	system("gnuplot PCAORIG.gp");
end	
