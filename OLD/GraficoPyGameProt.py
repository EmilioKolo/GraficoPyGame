###Libraries

import pygame
import os

##############################FUNCIONES##############################

def abrirTXT(nombre,extension='csv'):
    F = open(os.path.dirname(os.path.abspath(__file__)) + '\\' + str(nombre) + '.' + str(extension),'r');
    L = F.readlines();
    F.close();
    return L

def ordenarL(L,orden):
    ordL = [0 for i in orden];
    for i in range(len(orden)):
        ordL[i] = L[orden[i]];
    return ordL

def sacarFeno(nom,num,ext='csv',headers=True,sep=';'):
    '''
    num son numeros del 0-9 para referirse a
    0-TODO en la base de datos
    1-TODO en poblacion o papers
    2-Patogenico
    3-Poblacional o Polimorfismo
    4-Indefinido
    5-
    6-
    7-
    8-
    9-
    '''
    L = abrirTXT(nom,extension=ext);
    head = L[0].rstrip('\n').split(sep);
    
    Feno = [];
    for i in L[int(headers):]:
        l = i.rstrip('\n').split(sep);
        if str(l[num+2]) == str(head[num+2]):
            Feno.append(l[0]);
    return Feno

def sacarTodo(nom,ext='csv',headers=True,sep=';'):
    L = abrirTXT(nom,extension=ext);
    Todo = [];
    for i in L[int(headers):]:
        l = i.rstrip('\n').split(sep);
        Todo.append(l[0]);
    return Todo

##############################VARIABLES##############################

#Colores predefinidos
blanco = (255,255,255);
negro = (0,0,0);
rojo = (200,0,0);
verde = (0,255,0);
azul = (0,0,255);
amarillo = (200,200,0);
verdeazul = (0,255,255);
violeta = (255,0,255);
naranja = (255,125,0);
gris = (150,150,150);
verdeoscuro = (0,125,0);
grisoscuro = (50,50,50);
rojooscuro = (100,0,0);
violetaoscuro = (125,0,125);
azuloscuro = (0,0,125);
amarillooscuro = (100,100,0);

#Colores como variables, para poner en funciones
colores = [amarillooscuro,verdeoscuro,rojo,grisoscuro,azuloscuro,grisoscuro,grisoscuro,grisoscuro,grisoscuro];
colorUTRintron = (175,175,175)

#Margen que se deja
margen = [100,100,100,100];

#Nombre del archivo de datos
nomDatos = 'DatosPyGameProt';

#Lista de columnas para extraer del archivo de datos
listaCol = [4,3,2,1,0];

#Ancho de las lineas de GV y de los exones
ancho = 96;
#Espaciado entre series de lineas
espaciado = 32;

#Lista de filas con representacion de la proteina
L_T = [6,7,8,9,10,11];

#############################INFO NKX2.5#############################

largoNKX = 975;
colorexones = negro;
w = 1;
largolinea = 10;
largoE1 = 335;
largoE2 = 640;
#Posiciones de inicio/fin de dominios en cDNA
Tinman = (28,63);
Homeodominio = (412,591);
NK2especifico = (634,702);
TyrRich = (694,786);
Helice1 = (436,474);
Helice2 = (493,528);
Helice3 = (535,582);
NLS1 = (406,429);
NLS2 = (574,597);
motivoSUMO = (151,162);

#Exones (version vieja)
exon1 = (colorexones,margen[0]+largoE2,margen[0],largoE1,ancho);
exon2 = (colorexones,margen[0],margen[0],largoE2,ancho);

#Proteina
prot = (gris,margen[0],margen[0],largoNKX,ancho/2)

#Dominios confirmados
TN = (amarillo,margen[0]+largoNKX-Tinman[1],margen[0],Tinman[1]-Tinman[0]+1,ancho/2);
HD = (rojooscuro,margen[0]+largoNKX-Homeodominio[1],margen[0],Homeodominio[1]-Homeodominio[0]+1,ancho/2);
NK2SD = (violetaoscuro,margen[0]+largoNKX-NK2especifico[1],margen[0],NK2especifico[1]-NK2especifico[0]+1,ancho/2);

#Dominios secundarios
YRD = (azuloscuro,margen[0]+largoNKX-TyrRich[1],margen[0]+ancho/4,TyrRich[1]-TyrRich[0]+1,ancho/4);
H1 = (rojo,margen[0]+largoNKX-Helice1[1],margen[0]+ancho/4,Helice1[1]-Helice1[0]+1,ancho/4);
H2 = (rojo,margen[0]+largoNKX-Helice2[1],margen[0]+ancho/4,Helice2[1]-Helice2[0]+1,ancho/4);
H3 = (rojo,margen[0]+largoNKX-Helice3[1],margen[0]+ancho/4,Helice3[1]-Helice3[0]+1,ancho/4);
NLS1 = (azuloscuro,margen[0]+largoNKX-NLS1[1],margen[0]+ancho/2,NLS1[1]-NLS1[0]+1,ancho/4);
NLS2 = (azuloscuro,margen[0]+largoNKX-NLS2[1],margen[0]+ancho/2,NLS2[1]-NLS2[0]+1,ancho/4);

nkxdraw = [exon2,exon1];
dominios1 = [TN,HD,NK2SD];
dominios2 = [YRD,H1,H2,H3,NLS1,NLS2];

##############################PANTALLA###############################

#Saco la lista de variantes genicas para graficar
ListaFeno = ordenarL([sacarFeno(nomDatos,i) for i in range(len(listaCol))],listaCol) + [sacarTodo(nomDatos)];

#Arranco los graficos
pygame.init();

#Tamaño de pantalla
Display = pygame.display.set_mode((1200,2500));
Display.fill(blanco);

### ACA VA LA PARTE QUE SE IMPRIME ###

#Por cada valor de la lista L_T
for t in L_T:
    #Dibujo los exones
    #pygame.draw.rect(Display,exon1[0],(exon1[1],exon1[2]+(t*ancho)+(t*espaciado),exon1[3],exon1[4]),w);
    #pygame.draw.rect(Display,exon2[0],(exon2[1],exon2[2]+(t*ancho)+(t*espaciado),exon2[3],exon2[4]),w);

    #Dibujo la proteina y marco el corte entre exones
    pygame.draw.rect(Display,prot[0],(prot[1],prot[2]+(t*ancho)+(t*espaciado),prot[3],prot[4]));
    pygame.draw.line(Display,prot[0],
                     (margen[0]+largoE2,margen[0]+(t*ancho)+(t*espaciado)-largolinea),
                     (margen[0]+largoE2,margen[0]+(t*ancho)+(t*espaciado)),w);

#Agrego las variantes
for i in range(2*len(ListaFeno)):
    for j in ListaFeno[i%len(ListaFeno)]:
        pygame.draw.line(Display, colores[i%len(ListaFeno)],
                        (margen[0]+largoNKX-int(j),margen[0]+(i*ancho)+(i*espaciado)),
                        (margen[0]+largoNKX-int(j),margen[0]+(i*ancho)+(i*espaciado)+ancho/2-1), 1);


#Dibujo las proteinas del final con dominios
#Primero con los dominios principales
pygame.draw.rect(Display,prot[0],(prot[1],prot[2]+(12*ancho)+(12*espaciado),prot[3],prot[4]));
pygame.draw.line(Display,prot[0],
                (margen[0]+largoE2,margen[0]+(12*ancho)+(12*espaciado)-largolinea),
                (margen[0]+largoE2,margen[0]+(12*ancho)+(12*espaciado)),w);
for i in dominios1:
    pygame.draw.rect(Display,i[0],(i[1],i[2]+(12*ancho)+(12*espaciado),i[3],i[4]));

#Despues con dominios menos importantes tambien
pygame.draw.rect(Display,prot[0],(prot[1],prot[2]+(13*ancho)+(13*espaciado),prot[3],prot[4]));
pygame.draw.line(Display,prot[0],
                (margen[0]+largoE2,margen[0]+(13*ancho)+(13*espaciado)-largolinea),
                (margen[0]+largoE2,margen[0]+(13*ancho)+(13*espaciado)),w);
for i in dominios1:
    pygame.draw.rect(Display,i[0],(i[1],i[2]+(13*ancho)+(13*espaciado),i[3],i[4]));
for i in dominios2:
    pygame.draw.rect(Display,i[0],(i[1],i[2]+(13*ancho)+(13*espaciado),i[3],i[4]));

#Por ultimo solo la forma de la proteina
pygame.draw.rect(Display,prot[0],(prot[1],prot[2]+(14*ancho)+(14*espaciado),prot[3],prot[4]));
pygame.draw.line(Display,prot[0],
                (margen[0]+largoE2,margen[0]+(14*ancho)+(14*espaciado)-largolinea),
                (margen[0]+largoE2,margen[0]+(14*ancho)+(14*espaciado)),w);

### ACA VA LA PARTE QUE SE IMPRIME ###

#Esto es para guardarlo y cerrar
pygame.display.update();
Display.set_alpha(255);
pygame.image.save(Display, 'FiguraProt.png');
pygame.quit()
