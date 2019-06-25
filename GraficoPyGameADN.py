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
verdeoscuro = (0,125,0);
grisoscuro = (50,50,50);
amarillooscuro = (100,100,0);

colores = [amarillooscuro,verdeoscuro,rojo,grisoscuro,azul,grisoscuro,grisoscuro,grisoscuro,grisoscuro];
colorUTRintron = (175,175,175)
#Colores como variables, para poner en funciones

margen = [100,100,100,100];
#Margen que se deja

nomDatos = 'DatosPyGame';
#Nombre del archivo de datos

ListaFeno = ordenarL([sacarFeno(nomDatos,i) for i in range(5)],[4,3,2,1,0]) + [sacarTodo(nomDatos)];
#Lista de variantes genicas para graficar

ancho = 96;
#Ancho de las lineas
espaciado = 64;
#Espaciado entre datos

#############################INFO NKX2.5#############################

posy = margen[0]; #  + (x*ancho) + (x*espaciado)
dif = 173231702;
largoNKX = 3208;
largoISO4 = 1813;
difISO4 = 1405;
colorexones = negro;
w = 10;

utr = (colorUTRintron,
       margen[0]+173232104-dif,
       posy+(ancho/3),
       largoNKX,
       ancho/3,
       0);
exon2 = (colorexones,
         margen[0]+173232569-dif,
         posy,
         641,
         ancho,
         w);
exon1 = (colorexones,
         margen[0]+173234750-dif,
         posy,
         334,
         ancho,
         w);

#UTR 3' del exon 2
UTR3 = (colorexones,
        margen[0]+173232104-dif,
        posy,
        465,
        ancho,
        w);

#UTR 5' del exon 1
UTR5 = (colorexones,
        margen[0]+173235084-dif,
        posy,
        229,
        ancho,
        w);

#Variantes de isoforma 2
exon2ISO2 = (colorexones,
            margen[0]+173233343-dif,
            posy,
            122,
            ancho,
            w);
UTR3ISO2 = (colorexones,
            margen[0]+173232104-dif,
            posy,
            1239,
            ancho,
            w);

#Variantes de isoforma 3
exon2ISO3 = (colorexones,
            margen[0]+173233497-dif,
            posy,
            5,
            ancho,
            w);
UTR3ISO3 = (colorexones,
            margen[0]+173232104-dif,
            posy,
            1393,
            ancho,
            w);

#Variantes de isoforma 4 (predicha)
utrISO4 = (colorUTRintron,
           margen[0]+173233509-dif,
           posy+(ancho/3),
           largoISO4,
           ancho/3,
           0);
UTR5ISO4 = (colorexones,
            margen[0]+173235084-dif,
            posy,
            238,
            ancho,
            w);
UTR3ISO4 = (colorexones,
            margen[0]+173233509-dif,
            posy,
            537,
            ancho,
            w);

exon2ISO4 = (colorexones,
            margen[0]+173234046-dif,
            posy,
            95,
            ancho,
            w);

#Listas para dibujar cada isoforma
nkxdraw = [utr, UTR3, UTR5, exon2, exon1];

nkxISO2 = [utr, UTR3ISO2, UTR5, exon2ISO2, exon1];
nkxISO3 = [utr, UTR3ISO3, UTR5, exon2ISO3, exon1];

nkxISO4 = [utrISO4, UTR3ISO4, UTR5ISO4, exon2ISO4, exon1];

##############################PANTALLA###############################

pygame.init();

Display = pygame.display.set_mode((4000,2650));
Display.fill(blanco);

#Primero dibujo el gen al fondo
t = 2*len(ListaFeno);

for n in nkxdraw:
    if bool(n[5]):
        pygame.draw.rect(Display,blanco,(n[1],n[2]+(t*ancho)+(t*espaciado),n[3],n[4]));
    pygame.draw.rect(Display,n[0],(n[1],n[2]+(t*ancho)+(t*espaciado),n[3],n[4]),n[5]);

#Despues pongo todas las variantes
for i in range(2*len(ListaFeno)):
    b = i >= len(ListaFeno);
    if b:
        #Primero dibujo el esquema de nkx2-5
        for n in nkxdraw:
            if bool(n[5]):
                pygame.draw.rect(Display,blanco,(n[1],n[2]+(i*ancho)+(i*espaciado),n[3],n[4]));
            pygame.draw.rect(Display,n[0],(n[1],n[2]+(i*ancho)+(i*espaciado),n[3],n[4]),n[5]);
        #Despues agrego las variantes
        for j in ListaFeno[i%len(ListaFeno)]:
            pygame.draw.line(Display, colores[i%len(ListaFeno)],
                             (margen[3]+int(j)-dif,margen[0]+(i*ancho)+(i*espaciado)),
                             (margen[3]+int(j)-dif,margen[0]+(i*ancho)+ancho+(i*espaciado)), 1);
    else:
        #Aca solo agrego las variantes
        for j in ListaFeno[i]:
            pygame.draw.line(Display, colores[i],
                             (margen[3]+int(j)-dif,margen[0]+(i*ancho)+(i*espaciado)),
                             (margen[3]+int(j)-dif,margen[0]+(i*ancho)+ancho+(i*espaciado)), 1);

#Despues agrego las isoformas
t = t + 1;
for n in nkxISO2:
    if bool(n[5]):
        pygame.draw.rect(Display,blanco,(n[1],n[2]+(t*ancho)+(t*espaciado),n[3],n[4]));
    pygame.draw.rect(Display,n[0],(n[1],n[2]+(t*ancho)+(t*espaciado),n[3],n[4]),n[5]);

t = t + 1;
for n in nkxISO3:
    if bool(n[5]):
        pygame.draw.rect(Display,blanco,(n[1],n[2]+(t*ancho)+(t*espaciado),n[3],n[4]));
    pygame.draw.rect(Display,n[0],(n[1],n[2]+(t*ancho)+(t*espaciado),n[3],n[4]),n[5]);

t = t + 1;
for n in nkxISO4:
    if bool(n[5]):
        pygame.draw.rect(Display,blanco,(n[1],n[2]+(t*ancho)+(t*espaciado),n[3],n[4]));
    pygame.draw.rect(Display,n[0],(n[1],n[2]+(t*ancho)+(t*espaciado),n[3],n[4]),n[5]);

pygame.display.update();
Display.set_alpha(255);
pygame.image.save(Display, 'FiguraADN.png');
pygame.quit()
