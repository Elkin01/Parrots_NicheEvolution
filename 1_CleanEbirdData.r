rm(list=ls())

library(raster)
library(sf)
library(terra)

alt=raster("~/OneDrive - City University of New York/Documentos/IAvH Humboldt/Montañas/altitude.tif")

alt2=aggregate(alt, 10, FUN="mean")

mundo=st_read('~/Library/CloudStorage/OneDrive-CityUniversityofNewYork/Documentos/Estadistica R SIG/Capas SIG/Shapes/Political_World_Borders/TM_WORLD_BORDERS_SIMPL-0.3/TM_WORLD_BORDERS_SIMPL-0.3.shp', 'TM_WORLD_BORDERS_SIMPL-0.3')
mundo <- mundo[1] 
mundo=as_Spatial(mundo)

birdlife=read.delim('~/OneDrive - City University of New York/CUNY PhD/Thesis/Cap 1 Macro Montanas/Tablas, filogenias y bases de datos/Final_Birdfile_V5_2021.txt', header=T)

registros=list.files("~/Dropbox (AMNH)/parrot_morpho_evolution/ebird_data_with_headers", pattern=".txt")

sinonimias=read.delim('~/Dropbox (AMNH)/parrot_morpho_evolution/ebird_data_with_headers/Filtered Ebird records/SPMatchParrotsMaps.txt', header=T)

setwd('~/Dropbox (AMNH)/parrot_morpho_evolution/ebird_data_with_headers')

tablaResumen=data.frame()

pdf('eBird_filter.pdf', width=25, height=10)

for (i in 333:length(registros))

{

  tabla=read.delim(registros[i], header=T)

  Species=sinonimias[match(unique(tabla$SCIENTIFIC.NAME), sinonimias[,'Species']),2]

  if(!is.na(Species))
  
    {

      Species_map=st_read("~/OneDrive - City University of New York/Documentos/Estadistica R SIG/Capas SIG/Shapes/Mapas BirdLife 2020/Shapefiles BirdLife 2020", strsplit(as.character(Species), split="[.]")[[1]][1]) # Just reading the map
     
      mapa <- Species_map[Species_map$ORIGIN == 1 | Species_map$ORIGIN == 2,] ## Leaving just distrbutions considered as "Native" and "Reintroduced " according to birdlife
      
      mapa=as_Spatial(mapa) #changing format to "SpatialPolygonsDataFrame"

      #mapa=st_combine(mapa)
      
      #mapa=st_buffer(mapa, 100000)
      
      #mapa<-st_make_valid(mapa)
      
      #sf_use_s2(FALSE)
      
      mapa1=buffer(mapa, 100000) #Generating a buffer of 1 degree
      
      coordinates(tabla)=c('LONGITUDE', 'LATITUDE')
 
      mapRaster=crop(alt2, mapa1@bbox[c(1, 3, 2, 4)])
 
      mapRaster=mask(mapRaster, mapa1)
      mapRaster[which(!is.na(values(mapRaster)))]<-1
 
      extraccion=extract(mapRaster, tabla)
 
      tabla=data.frame(tabla, extraccion)
 
      tabla2=tabla[which(!is.na(tabla$extraccion)),]
 
 
      if (nrow(tabla2)<5)
        
        {
                
        tabla2=as.data.frame(matrix(ncol=ncol(tabla2), nrow=100))
        colnames(tabla2)=colnames(tabla)
               
        randomCoor=spsample(mapa1, 100, type="random")
        
        tabla2=data.frame(tabla2, xyFromCell(mapRaster, extract(mapRaster, randomCoor, cellnumbers=TRUE)[,1]), Celda=extract(mapRaster, randomCoor, cellnumbers=TRUE)[,1])
        
        tabla2[,'LONGITUDE']<-tabla2$x
        tabla2[,'LATITUDE']<-tabla2$y
        
        
        NumRegistros=as.data.frame(tapply(tabla2$extraccion, INDEX = list(tabla2$Celda), 'sum'))
        colnames(NumRegistros)=c('NumRegistros')
        NumRegistros=data.frame(Celda=rownames(NumRegistros),NumRegistros)
        
        listaCeldas=as.numeric(unique(NumRegistros$Celda))
        
        }
 
      if (nrow(tabla2)>4)
    
       {
 
        tabla3=tabla2

        coordinates(tabla3)=c('LONGITUDE', 'LATITUDE')
 
        tabla2=data.frame(tabla2, xyFromCell(mapRaster, extract(mapRaster, tabla3, cellnumbers=TRUE)[,1]), Celda=extract(mapRaster, tabla3, cellnumbers=TRUE)[,1])
 
        NumRegistros=as.data.frame(tapply(tabla2$extraccion, INDEX = list(tabla2$Celda), 'sum'))
        colnames(NumRegistros)=c('NumRegistros')
        NumRegistros=data.frame(Celda=rownames(NumRegistros),NumRegistros)
 
        listaCeldas=as.numeric(unique(NumRegistros$Celda))
 
       }
 
    }


  if(is.na(Species))
  
    {
  
     if(nrow(tabla)<5)
      {
       tablaResumen[i,1]<-unique(tabla$SCIENTIFIC.NAME)
       tablaResumen[i,2]<-Species
       tablaResumen[i,3]<-registros[i]
       tablaResumen[i,4]<-nrow(tabla)
       tablaResumen[i,5]<-0
      } else {next}
  
     if(nrow(tabla)>4)
     
      {
  
       tabla2=tabla
  
       tabla3=tabla2
  
       coordinates(tabla3)=c('LONGITUDE', 'LATITUDE')
  
       mapRaster=crop(alt2, c(min(tabla$LONGITUDE), max(tabla$LONGITUDE), min(tabla$LATITUDE), max(tabla$LATITUDE)))
  
      tabla2=data.frame(tabla2, xyFromCell(mapRaster, extract(mapRaster, tabla3, cellnumbers=TRUE)[,1]), Celda=extract(mapRaster, tabla3, cellnumbers=TRUE)[,1])
  
      tabla2[,'extraccion']<-1 #Se adiciona un uno en esta columna para hacer el conteo de filas como el número de registros en el siguiente comando. A esta columna se le adiciona el nombre de "extraccion" para hacer esta linea de comando equivalente al "if" alternativo cuando la especie tiene mapa. 
  
      NumRegistros=as.data.frame(tapply(tabla2$extraccion, INDEX = list(tabla2$Celda), 'sum'))
      colnames(NumRegistros)=c('NumRegistros')
      NumRegistros=data.frame(Celda=rownames(NumRegistros),NumRegistros)
  
      listaCeldas=as.numeric(unique(NumRegistros$Celda))
     }
  
    }


  tablaF=data.frame()
 
  for (j in 1:length(listaCeldas))
   
  {
   
   tabla4=subset(tabla2, Celda==listaCeldas[j])
   
   tabla4=tabla4[!duplicated(tabla4$OBSERVER.ID), ]
   
    if (dim(tabla4)[1]<2) {next}
   
    if (dim(tabla4)[1]>1)
     
      {
       tablaF=rbind(tablaF, tabla4)
      }
   

  }


###################################################
### depuracion 2, solo un registro por coordenada de  1km
###################################################

  tabla3=tablaF[,c('LONGITUDE', 'LATITUDE')]
  coordinates(tabla3)=c('LONGITUDE', 'LATITUDE')
  
  tablaF=data.frame(tablaF, xyFromCell(alt, extract(alt, tabla3, cellnumbers=TRUE)[,1]), Celda2=extract(alt, tabla3, cellnumbers=TRUE)[,1])
  
  tablaF=tablaF[!duplicated(tablaF[,c('SCIENTIFIC.NAME', 'x.1', 'y.1')]),]

###################################################


  tablaResumen[i,1]<-unique(tabla$SCIENTIFIC.NAME)
  tablaResumen[i,2]<-Species
  tablaResumen[i,3]<-registros[i]
  tablaResumen[i,4]<-nrow(tabla)
  tablaResumen[i,5]<-nrow(tablaF)
 
  write.table(tablaF, file=paste(1, registros[i], sep = '_'), sep = '\t')
 
   ###archivos temporales solo para hacer plots de chequeo
   
 
   par(mfrow=c(1,3))
   tablaF1=tablaF
   coordinates(tablaF1)=c('LONGITUDE', 'LATITUDE')
   plot(mundo, main=Species)
   plot(tablaF1, pch=21, col='red', add=T)
   plot(mapRaster)
   
  
   plot(tablaF1, add=T)
   
   
   NumRegistros=as.data.frame(tapply(tablaF$extraccion, INDEX = list(tablaF$Celda), 'sum'))
   colnames(NumRegistros)=c('NumRegistros')
   NumRegistros=data.frame(Celda=rownames(NumRegistros),NumRegistros)
   
   mapRaster[]<-NA
   mapRaster[as.numeric(NumRegistros$Celda)]<-NumRegistros$NumRegistros
   plot(mapRaster)
   plot(mapa, add=T)
 
}

dev.off()

write.table(tablaResumen, file='SummaryEbirdFilteredRecords.txt', sep = '\t')

#tabla=tabla[-which(is.na(tabla$LATITUDE) | is.na(tabla$LONGITUDE)),]
length(which(duplicated(tabla2[,c('LATITUDE', 'LONGITUDE', 'OBSERVER.ID')])==TRUE))





