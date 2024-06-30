
# Generales
import os 
import time 
import copy 
import math 
import sys 
import pandas as pd 
import numpy as np 

# Para display
from random import shuffle 

# Graficos
import matplotlib.pyplot as plt 
import seaborn as sns 
import pandas as pd
import pygame


'''
Funciones para generar un grafico de visualizacion de variantes puntuales
'''

#################################### VARIABLES ####################################

### Variables main()

# Path main depende de posicion actual
path_git_main = os.path.dirname(os.path.abspath(__file__)); 
path_output_dump_main = 'D:\\ArchivosDoctorado\\Output_dump\\'; 
path_in_main = path_git_main + ''; ### COMPLETAR
path_out_main = path_output_dump_main + ''; ### COMPLETAR

nom_archivo_variantes = 'tabla_variantes'; 
path_archivo_variantes = path_git_main; 

# Columnas seleccionadas, las primeras 3 son identificadores, 
# la 4ta es para seleccionar transition/transversion, la 5ta es para seleccionar missense
# la 6ta es para seleccionar pathogenic y healthy, la 7ma es para seleccionar frecuencias
nom_columnas = ["DNA", "cDNA", "Protein", "GV Type", "Protein effect", "Classification", "Frequency per 1000 ind"]; 

nom_out_main = 'FiguraDistribucion'; 


#################################### FUNCIONES ####################################


def _main():
    # Funcion para probar funciones en ejecucion del archivo

    pipeline_grafico(nom_archivo_variantes, nom_out_main, nom_columnas, path_gvs=path_archivo_variantes, path_out=path_archivo_variantes); 

    return ''


def pipeline_grafico(nom_gvs, nom_out, nom_cols, path_gvs='', path_out=''): 
    '''Genera el grafico de GVs por ventana desde el archivo de variantes y las columnas correspondientes'''

    # Abro el archivo y guardo l_head
    m_gvs, l_head = abrir_csv(nom_gvs, path_arch=path_gvs, con_headers=True, devolver_header=True); 

    # Uso la funcion seleccionar_columnas_filtro para seleccionar los datos correctos y transformarlos
    m_gvs_adn, m_gvs_prot = seleccionar_columnas_filtro(m_gvs, l_head, nom_cols); 

    # Uso la funcion generar_grafico para hacer graficos con los datos filtrados
    generar_grafico_adn(m_gvs_adn, nom_out + '_adn', path_out=path_out); 
    generar_grafico_prot(m_gvs_prot, nom_out + '_prot', path_out=path_out); 
    
    return m_gvs_adn, m_gvs_prot 


### Funciones principales del pipeline


def generar_grafico_adn(m_gvs, nom_out, path_out=''):
    '''Funcion que hace un grafico en base a los datos ya procesados'''

    ancho = 96;
    # Ancho de las lineas
    espaciado = 64;
    # Espaciado entre datos

    # Colores predefinidos
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
    # Colores como variables, para poner en funciones

    margen = [100,100,100,100];

    ### INFO NKX2.5

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

    ### PANTALLA
    pygame.init();
    Display = pygame.display.set_mode((4000,2650));
    Display.fill(blanco);
    #Primero dibujo el gen al fondo
    t = 2*len(m_gvs);

    for n in nkxdraw:
        if bool(n[5]):
            pygame.draw.rect(Display,blanco,(n[1],n[2]+(t*ancho)+(t*espaciado),n[3],n[4]));
        pygame.draw.rect(Display,n[0],(n[1],n[2]+(t*ancho)+(t*espaciado),n[3],n[4]),n[5]);

    #Despues pongo todas las variantes
    for i in range(2*len(m_gvs)):
        b = i >= len(m_gvs);
        if b:
            #Primero dibujo el esquema de nkx2-5
            for n in nkxdraw:
                if bool(n[5]):
                    pygame.draw.rect(Display,blanco,(n[1],n[2]+(i*ancho)+(i*espaciado),n[3],n[4]));
                pygame.draw.rect(Display,n[0],(n[1],n[2]+(i*ancho)+(i*espaciado),n[3],n[4]),n[5]);
            #Despues agrego las variantes
            for j in m_gvs[i%len(m_gvs)]:
                pygame.draw.line(Display, colores[i%len(m_gvs)],
                                (margen[3]+int(j)-dif,margen[0]+(i*ancho)+(i*espaciado)),
                                (margen[3]+int(j)-dif,margen[0]+(i*ancho)+ancho+(i*espaciado)), 1);
        else:
            #Aca solo agrego las variantes
            for j in m_gvs[i]:
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
    pygame.image.save(Display, os.path.join(path_out, nom_out+'.png'));
    pygame.quit()


def generar_grafico_prot(m_gvs, nom_out, path_out=''):
    '''Funcion que hace un grafico en base a los datos ya procesados'''

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
    #Lista de columnas para extraer del archivo de datos
    listaCol = [4,3,2,1,0];

    #Ancho de las lineas de GV y de los exones
    ancho = 96;
    #Espaciado entre series de lineas
    espaciado = 32;

    #Lista de filas con representacion de la proteina
    L_T = [6,7,8,9,10,11];

    ### INFO NKX2.5

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

    # Arranco los graficos
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
    for i in range(2*len(m_gvs)):
        for j in m_gvs[i%len(m_gvs)]:
            pygame.draw.line(Display, colores[i%len(m_gvs)],
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
    pygame.image.save(Display, os.path.join(path_out, nom_out+'.png'));
    pygame.quit()


def seleccionar_columnas_filtro(m_gvs, l_head, nom_cols): 
    '''Selecciona columnas en nom_cols y filtra y procesa las cosas'''

    # Inicializo las matrices que se devuelven
    m_out_adn = [[],[],[],[],[],[]]; 
    m_out_prot = [[],[],[],[],[],[]]; 
    # Defino ids de nom_cols
    l_ids = []; 
    for k in range(len(nom_cols)): 
        # Agrego el id correspondiente a l_ids
        l_ids.append(l_head.index(nom_cols[k])); 
    # Recorro m_gvs
    for i in range(len(m_gvs)): 
        curr_gv = m_gvs[i]; 
        # Selecciono transition/transversion
        if curr_gv[l_ids[3]].lower() in ['transition', 'transversion']: 
            # Transformo el identificador de ADN en posicion
            pos_adn = procesar_adn(curr_gv[l_ids[0]]); 
            # Agrego pos_adn a m_out_adn[4] y m_out_adn[5]
            m_out_adn[4].append(int(pos_adn)); 
            m_out_adn[5].append(int(pos_adn)); 
            # Selecciono pathogenic
            if curr_gv[l_ids[5]].lower()=='pathogenic': 
                m_out_adn[2].append(int(pos_adn)); 
            # Selecciono benignas y frecuencias altas
            if (curr_gv[l_ids[5]].lower() in ['likely benign', 'benign']) or (curr_gv[l_ids[6]]!='' and float(curr_gv[l_ids[6]])>=0.1): 
                m_out_adn[1].append(int(pos_adn)); 
            # Selecciono variantes en exones 1 o 2 por su efecto en la proteina
            if curr_gv[l_ids[4]].lower() in ['missense', 'nonsense', 'stop loss', 'synonymous']: 
                # Transformo el identificador de ADNc en posicion de proteina
                pos_adnc = procesar_adn(curr_gv[l_ids[1]]); 
                pos_prot = int(pos_adnc); 
                m_out_prot[4].append(int(pos_prot)); 
                m_out_prot[5].append(int(pos_prot)); 
                # Selecciono pathogenic
                if curr_gv[l_ids[5]].lower()=='pathogenic': 
                    m_out_prot[2].append(int(pos_prot)); 
                # Selecciono benignas y frecuencias altas
                if (curr_gv[l_ids[5]].lower() in ['likely benign', 'benign']) or (curr_gv[l_ids[6]]!='' and float(curr_gv[l_ids[6]])>=0.1): 
                    m_out_prot[1].append(int(pos_prot)); 
    return m_out_adn, m_out_prot


### Funciones simples de procesamiento


def procesar_adn(id_adnc):
    '''Recibe un identificador de ADNc en formato HGVS y devuelve solo la posicion'''

    # Inicializo el texto que se devuelve
    str_out = ''; 
    # Borro el c. del principio y el cambio N>N del final
    str_out = id_adnc[2:-3]; 
    return int(str_out)


def procesar_prot(id_prot):
    '''Recibe un identificador de proteina en formato HGVS y devuelve solo la posicion'''

    # Inicializo el texto que se devuelve
    str_out = ''; 
    # Veo si id_prot tiene parentesis
    if "(" in id_prot:
        # Borro el principio y el final
        str_out = id_prot[6:-4]; 
    else:
        # Borro el principio y el final
        str_out = id_prot[5:-3]; 
    return int(str_out)


### Funciones para abrir y guardar archivos


def abrir_csv(nom_arch, path_arch='', sep=';', ext='.csv', con_headers=True, devolver_header=False):
    '''Abre archivos .csv y devuelve una lista de listas con las filas.
    con_headers define si hay header (se ignora por defecto).
    devolver_header permite devolver el header como segundo output'''

    # Inicializo la matriz que se devuelve
    M_out = []; 
    # Inicializo booleano y lista de headers
    header = con_headers; 
    l_head = []; 
    # Defino la direccion del archivo con nom_arch y path_arch
    if path_arch=='':
        filepath = nom_arch + ext; 
    else:
        filepath = os.path.join(path_arch, nom_arch + ext); 
    # Abro el archivo filepath
    with open(filepath, 'r') as f_out:
        print('Archivo ' + nom_arch + ' abierto.')
        # Recorro cada fila
        for curr_line in f_out.readlines():
            # Veo si tengo que pasar el header
            if header:
                # Transformo la linea en lista
                l_head = curr_line.rstrip().split(sep=sep); 
                # Paso header a False para no ver la proxima linea
                header = False; 
            else:
                # Transformo la linea en lista
                l_line = curr_line.rstrip().split(sep=sep); 
                # Cargo l_line en M_out
                M_out.append(l_line[:]); 
    # Si devolver_header es True, se devuelve M_out y l_head
    if devolver_header:
        return M_out, l_head
    return M_out


#################################### RESULTADOS ###################################

output_dump = []; 

if __name__=='__main__':
    output_dump.append(_main()); 

