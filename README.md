17 de Junio de 2013 - Jorge Velásquez-Tibatá - Instituto Alexander von Humboldt

###mxParallel

La función mxParallel tiene el propósito de servir como envoltura para la ejecución en paralelo de Maxent (Phillips et al. 2006) desde R, así como de guiar decisiones claves de modelamiento (e.g. selección de área de estudio, método de extrapolación, features etc.) por medio de cuadros de diálogo.

La ejecución comienza determinando como workspace el directorio donde residen las funciones mxParallel.R y sdm-functions.R y cargando mxParallel:

    setwd("C:/Users/RScripts")
    source("mxParallel.R")

Aunque el código cargara y de ser necesario instalará los paquetes necesarios para su ejecución, en el caso particular del paquete dismo deberá descargar el archivo .jar de Maxent manualmente y copiarlo en la carpeta java del paquete dismo la cual esta ubicada en el directorio de instalación de R, o en la biblioteca de paquetes de R que aparece en Mis documentos.

Luego ejecute la función:

    res <- mxParallel()


El objeto res guardará los parámetros de ejecución escogidos interactivamente durante la sesión de mxParallel. A continuación escogerá la carpeta de archivos de salida donde serán guardados los modelos y la carpeta que alberga las variables continuas (en formato ESRI ASCII con extensión .asc).

![image1](https://lh3.googleusercontent.com/-2jBJx14qnr8/Ub9jkyvChuI/AAAAAAAAADs/MGn4jl7EVdA/w338-h384-no/Imagen1.png)
![image2](https://lh3.googleusercontent.com/-ula8yjLihkg/Ub9jk3dwrGI/AAAAAAAAAD0/vwH4se_Kl0U/w339-h385-no/Imagen2.png)

En el siguiente cuadro de diálogo podrá escoger las variables a usar en el modelamiento. Utilice la tecla Ctrl para seleccionar variables individuales o la tecla Shift para escoger un rango de variables.

![image3](https://lh6.googleusercontent.com/-TmmIBTfjpg8/Ub9jlE07OTI/AAAAAAAAAEA/0ZoZUtBuVNE/w188-h463-no/Imagen3.png)

En caso que sus capas ambientales no tengan información de proyección asociada, mxParallel asumirá que no están proyectadas y que tienen un datum WGS84. A continuación deberá ingresar la dirección de un archivo CSV con los datos de ocurrencia de las especies a ser modeladas. El archivo CSV deberá tener en su primera columna el nombre de la especie, y en la segunda y tercera longitud y latitud (o coordenadas X y Y en el mismo sistema de proyección de las capas ambientales), respectivamente. Campos adicionales serán ignorados.

![image4](https://lh3.googleusercontent.com/-bByyRf5P8G0/Ub9jlXdgIkI/AAAAAAAAAEM/qbmzbj3Lb0E/w641-h516-no/Imagen4.png)

El programa internamente eliminará de la base de datos las especies que tengan menos de 10 registros. Para cambiar el umbral de registros necesarios para elaborar un modelo, deberá hacerlo manualmente en la línea 33 del código de mxParallel. 

Ahora seleccione las especies que desea modelar:

![image5](https://lh6.googleusercontent.com/-HiUbnyelRjI/Ub9jlboYxbI/AAAAAAAAAEQ/r14bxGU8ZTY/w221-h766-no/Imagen5.png)

Los features permitidos de Maxent:

![image6](https://lh3.googleusercontent.com/-U9TPImNU_RQ/Ub9jlhN4GSI/AAAAAAAAAEU/x_Iogmvq9f4/w187-h208-no/Imagen6.png)

Las opciones de extrapolación a nuevos ambientes. En caso de escoger la opción none y otra, la opción none tendrá precedencia:

![image7](https://lh5.googleusercontent.com/-NArzHgvw7-M/Ub9jloznsPI/AAAAAAAAAEY/MTErNXhYWiU/w189-h168-no/Imagen7.png)

Los umbrales deseados para convertir los modelos de escala continua a binaria:

![image8](https://lh4.googleusercontent.com/-HIJ1Wbd1K08/Ub9jlxBjznI/AAAAAAAAAEk/NGPTJbUpx-w/w250-h190-no/Imagen8.png)

La opción "Custom percentile(s)" le permite escoger uno o varios umbrales basados en percentiles de probabilidad de ocurrencia de los datos de entrenamiento. Ingrese los valores deseados separados por comas:

![image9](https://lh6.googleusercontent.com/-uE4VDMF9KFU/Ub9jly0kkKI/AAAAAAAAAEo/m22WI_xOZHk/w336-h183-no/Imagen9.png)

A continuación seleccione si desea realizar evaluación de los modelos por medio de k-fold partitioning, y si desea cortar los modelos usando una regla de parche. La regla de parche lo que hace es retener como parte de un modelo binario (después de aplicar un umbral), aquellos parches (bloques contiguos de distribución predicha) en los que hay al menos una ocurrencia.

Luego seleccione el área de interés. Esta es el área de la cual se derivará el background para el modelo de Maxent. Existen tres opciones: Raster Extent, que toma el background de la extensión completa de las capas ambientales; Convex Hull, que toma el background del área definida por el polígono mínimo convexo que encierra las ocurrencias de cada especie y Regions, que toma el background de las regiones definidas por un shapefile (por ejemplo de provincias biogeográficas o ecoregiones) en las que hay datos de ocurrencia de cada especie. El área de interés en las opciones Convex Hull y Regions será distinta para cada especie, por lo que el tiempo de computación será en general mayor que en la opción Raster Extent.

![image10](https://lh6.googleusercontent.com/-98I0NiioWC0/Ub9jkngDYpI/AAAAAAAAADc/pejpDGJP3U0/w188-h164-no/Imagen10.png)

Finalmente podrá escoger si desea un background tomado aleatoriamente, o definido por un archivo con puntos que indican la intensidad de muestreo del grupo de interés (target background sensu  Phillips et al. 2009). 

![image11](https://lh4.googleusercontent.com/-BK2CBy2_Oqo/Ub9jkvCQHLI/AAAAAAAAADk/AnImsB4J7Gc/w188-h154-no/Imagen11.png)

Este código implementa por defecto la ejecución de Maxent utilizando 4 núcleos, los cuales son estándar en los procesadores Intel i5 e i7. Dependiendo del número de núcleos (físicos y virtuales) de su procesador podrá cambiar manualmente este número en la línea 120 del código. Los modelos serán guardados a medida que son realizados en la carpeta de salida especificada en formato grd, el cual se puede leer usando el freeware DIVA-GIS. Los modelos continuos son nombrados usando el nombre de la especie, los modelos binarios tienen un sufijo especificando el umbral utilizado (e.g. TP, MSS, ESS) y los modelos cortados usando la regla de parche tienen el sufijo adicional cut.

###Desarrollos futuros (Lista de tareas)
* Evaluar modelos usando jackniffe para especies con 5-10 datos de ocurrencia.
