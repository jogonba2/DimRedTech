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

## PCA Y LDA combinados (No son necesarios los 4 ultimos parametros), configurar a mano en la competicion ##
## Primero se aplica PCA proyectando hasta nº dimensiones de la imagen, después se aplica LDA sobre las 5 dimensiones que con PCA 
## den menor error ##
if(type=="PLC")
	minPCA     = 1;   # Ajustar dimension minima de PCA AQUI #
	maxPCA     = 256; # Ajustar dimension maxima de PCA AQUI #
	menoresErrorPCA = [];
	menoresDimPCA   = [];
	kk         = 1; # Configurar numero de vecinos del clasificador #
	while minPCA<maxPCA
		[avg,W] = pca(X,minPCA);
		proyX   = W'*X;
		proyY   = W'*Y;
		err     = knn(proyX,xl,proyY,yl,kk);
		if(minPCA<=5) # Coger las 5 dimensiones que den menor error con PCA #
		    menoresErrorPCA(end+1) = err;
		    menoresDimPCA(end+1) = minPCA;
		else
		    if(err<max(menoresErrorPCA))
		        i = find(menoresErrorPCA==max(menoresErrorPCA));
		        menoresErrorPCA(i) = err;
		        menoresDimPCA(i) = minPCA;
		    end
		end
		minPCA = minPCA + 1;
		#disp("Con error: "),disp(menoresErrorPCA),disp(" dimension: "),disp(menoresDimPCA)
	endwhile
	disp("Menores con PCA ya calculados: "),disp(menoresDimPCA),disp(menoresErrorPCA)
	# Calcular con los 5 menores de PCA el error con LDA y coger el menor #
	x = 1;
	j = 1;
	#menoresDimPCA = [10,20,30,40,50];
	menorError = 100;
	dimPCA = 0;
	dimLDA = 0;
	for x=1:5 # Con las 5 dimensiones de PCA que dan menor error #
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
	disp("Menor error con: "),disp(dimPCA),disp(dimLDA),disp(" : "),disp(menorError)
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
