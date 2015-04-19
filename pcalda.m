#!/usr/bin/octave -qf
if(nargin!=8)
    printf("Usage: pcalda.m <type> <trdata> <trlabels> <tsdata> <tslabels> <mink> <stepk> <maxk>\n");
    exit(1);
end

# Argumentos #
arg_list=argv();
type   = arg_list{1};
trdata = arg_list{2};
trlabs = arg_list{3};
tsdata = arg_list{4};
tslabs = arg_list{5};
mink   = str2num(arg_list{6});
stepk  = str2num(arg_list{7});
maxk   = str2num(arg_list{8});

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

fd = fopen("results.txt","w");
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

## PCA y LDA generando la gráfica ##
if(type=="PLA")
	fwrite(fd,"# This file has PCA and LDA results #\r\n");
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
