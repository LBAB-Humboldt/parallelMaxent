17 de Junio de 2013 - Jorge Vel�squez-Tibat� - Instituto Alexander von Humboldt

=mxParallel=

La funci�n mxParallel tiene el prop�sito de servir como envoltura para la ejecuci�n en paralelo de Maxent (Phillips et al. 2006) desde R, as� como de guiar decisiones claves de modelamiento (e.g. selecci�n de �rea de estudio, m�todo de extrapolaci�n, features etc.) por medio de cuadros de di�logo.

La ejecuci�n comienza determinando como workspace el directorio donde residen las funciones mxParallel.R y sdm-functions.R y cargando mxParallel:

{{{
setwd("C:/Users/RScripts")
source("mxParallel.R")
}}}

Aunque el c�digo cargara y de ser necesario instalar� los paquetes necesarios para su ejecuci�n, en el caso particular del paquete dismo deber� descargar el archivo .jar de Maxent manualmente y copiarlo en la carpeta java del paquete dismo la cual esta ubicada en el directorio de instalaci�n de R, o en la biblioteca de paquetes de R que aparece en Mis documentos.

Luego ejecute la funci�n:

{{{
res <- mxParallel()
}}}

El objeto res guardar� los par�metros de ejecuci�n escogidos interactivamente durante la sesi�n de mxParallel. A continuaci�n escoger� la carpeta de archivos de salida donde ser�n guardados los modelos y la carpeta que alberga las variables continuas (en formato ESRI ASCII con extensi�n .asc).

https://lh3.googleusercontent.com/-2jBJx14qnr8/Ub9jkyvChuI/AAAAAAAAADs/MGn4jl7EVdA/w338-h384-no/Imagen1.png
https://lh3.googleusercontent.com/-ula8yjLihkg/Ub9jk3dwrGI/AAAAAAAAAD0/vwH4se_Kl0U/w339-h385-no/Imagen2.png

En el siguiente cuadro de di�logo podr� escoger las variables a usar en el modelamiento. Utilice la tecla Ctrl para seleccionar variables individuales o la tecla Shift para escoger un rango de variables.

https://lh6.googleusercontent.com/-TmmIBTfjpg8/Ub9jlE07OTI/AAAAAAAAAEA/0ZoZUtBuVNE/w188-h463-no/Imagen3.png

En caso que sus capas ambientales no tengan informaci�n de proyecci�n asociada, mxParallel asumir� que no est�n proyectadas y que tienen un datum WGS84. A continuaci�n deber� ingresar la direcci�n de un archivo CSV con los datos de ocurrencia de las especies a ser modeladas. El archivo CSV deber� tener en su primera columna el nombre de la especie, y en la segunda y tercera longitud y latitud (o coordenadas X y Y en el mismo sistema de proyecci�n de las capas ambientales), respectivamente. Campos adicionales ser�n ignorados.

https://lh3.googleusercontent.com/-bByyRf5P8G0/Ub9jlXdgIkI/AAAAAAAAAEM/qbmzbj3Lb0E/w641-h516-no/Imagen4.png

El programa internamente eliminar� de la base de datos las especies que tengan menos de 10 registros. Para cambiar el umbral de registros necesarios para elaborar un modelo, deber� hacerlo manualmente en la l�nea 33 del c�digo de mxParallel. 

Ahora seleccione las especies que desea modelar:

https://lh6.googleusercontent.com/-HiUbnyelRjI/Ub9jlboYxbI/AAAAAAAAAEQ/r14bxGU8ZTY/w221-h766-no/Imagen5.png

Los features permitidos de Maxent:

https://lh3.googleusercontent.com/-U9TPImNU_RQ/Ub9jlhN4GSI/AAAAAAAAAEU/x_Iogmvq9f4/w187-h208-no/Imagen6.png

Las opciones de extrapolaci�n a nuevos ambientes. En caso de escoger la opci�n none y otra, la opci�n none tendr� precedencia:

https://lh5.googleusercontent.com/-NArzHgvw7-M/Ub9jloznsPI/AAAAAAAAAEY/MTErNXhYWiU/w189-h168-no/Imagen7.png

Los umbrales deseados para convertir los modelos de escala continua a binaria:

https://lh4.googleusercontent.com/-HIJ1Wbd1K08/Ub9jlxBjznI/AAAAAAAAAEk/NGPTJbUpx-w/w250-h190-no/Imagen8.png

La opci�n �Custom percentile(s)� le permite escoger uno o varios umbrales basados en percentiles de probabilidad de ocurrencia de los datos de entrenamiento. Ingrese los valores deseados separados por comas:

https://lh6.googleusercontent.com/-uE4VDMF9KFU/Ub9jly0kkKI/AAAAAAAAAEo/m22WI_xOZHk/w336-h183-no/Imagen9.png

A continuaci�n seleccione si desea realizar evaluaci�n de los modelos por medio de k-fold partitioning, y si desea cortar los modelos usando una regla de parche. La regla de parche lo que hace es retener como parte de un modelo binario (despu�s de aplicar un umbral), aquellos parches (bloques contiguos de distribuci�n predicha) en los que hay al menos una ocurrencia.

Luego seleccione el �rea de inter�s. Esta es el �rea de la cual se derivar� el background para el modelo de Maxent. Existen tres opciones: Raster Extent, que toma el background de la extensi�n completa de las capas ambientales; Convex Hull, que toma el background del �rea definida por el pol�gono m�nimo convexo que encierra las ocurrencias de cada especie y Regions, que toma el background de las regiones definidas por un shapefile (por ejemplo de provincias biogeogr�ficas o ecoregiones) en las que hay datos de ocurrencia de cada especie. El �rea de inter�s en las opciones Convex Hull y Regions ser� distinta para cada especie, por lo que el tiempo de computaci�n ser� en general mayor que en la opci�n Raster Extent.

https://lh6.googleusercontent.com/-98I0NiioWC0/Ub9jkngDYpI/AAAAAAAAADc/pejpDGJP3U0/w188-h164-no/Imagen10.png

Finalmente podr� escoger si desea un background tomado aleatoriamente, o definido por un archivo con puntos que indican la intensidad de muestreo del grupo de inter�s (target background sensu  Phillips et al. 2009). 

https://lh4.googleusercontent.com/-BK2CBy2_Oqo/Ub9jkvCQHLI/AAAAAAAAADk/AnImsB4J7Gc/w188-h154-no/Imagen11.png

Este c�digo implementa por defecto la ejecuci�n de Maxent utilizando 4 n�cleos, los cuales son est�ndar en los procesadores Intel i5 e i7. Dependiendo del n�mero de n�cleos (f�sicos y virtuales) de su procesador podr� cambiar manualmente este n�mero en la l�nea 120 del c�digo. Los modelos ser�n guardados a medida que son realizados en la carpeta de salida especificada en formato grd, el cual se puede leer usando el freeware DIVA-GIS. Los modelos continuos son nombrados usando el nombre de la especie, los modelos binarios tienen un sufijo especificando el umbral utilizado (e.g. TP, MSS, ESS) y los modelos cortados usando la regla de parche tienen el sufijo adicional cut.

==Desarrollos futuros (Lista de tareas)==
* Evaluar modelos usando jackniffe para especies con 5-10 datos de ocurrencia.