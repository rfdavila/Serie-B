from xml.etree.ElementTree import PI
import pandas as pd
from math import*
import matplotlib.pyplot as plt

class SerieB():
    def __init__(self,Pfac,Afac,z,re): #constructor para inicializar los valores de entrada
        
        self.Pfac=Pfac                 #Relacion paso/diametro(P/D)
        self.Afac=Afac                 #Coeficente de area (Ae/Ao)
        self.z=z                       #Numero de palas
        self.re=re                     #Numero de Reynolds
    
    def Write_Datos(self): #,runAllFunctions=True
        """
        if runAllFunctions:
            self.Curves_one()
        """
        
        
        
        Datos = [
            ['---PROPELLER---'],
            ['Curves Wageningen B-Screw'],
            ['P/D (0.5<P/D<1.4) :', self.Pfac, '[-]'],
            ['Ae/Ao (0.3<Ae/Ao<1.05) :', self.Afac, '[-]'],
            ['Z (number of Blades 2<z<7)', self.z, '[-]'],
            ['RE (number of Reynold 2E6<RE<1E9)', self.re, '[-]'],
        ]
        
        
        print(Datos)
        print(self.Curves_one())
        
        
        """"
        with open('Datos.txt','a+') as archivo:
            archivo.write("\n".join(Datos))
        """
        
        
        
       
    def Curves_one(self,Pfac,Afac,z,re):
        
        self.Pfac=Pfac                 
        self.Afac=Afac                 
        self.z=z                       
        self.re=re                     
        
        tabla=pd.DataFrame()                                          #Crea dataframe tabla donde se mostraran los resultados
        tabla['J']=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2]  #Crea columna en de tabla de los coef. de valores de avance
        tabla['KT']=0                                                 #Crea columna inicial=0 de valores de KT
        tabla['KQ']=0                                               #Crea columna inicial=0 de valores de 10KQ
        tabla['EFIC']=0
        dkt= pd.read_table('CKt.txt', sep='\s+')                      #Crea dataframe de valores importados de coeficiente de la regrecion para KT 
        dkq= pd.read_table('CKq.txt', sep='\s+')                      ##Crea dataframe de valores importados de coeficiente de la regrecion para 10KQ
        
        for i in range(len(tabla)):
            
            dkt['computkt']=(dkt.Tstuv)*((tabla.loc[i,'J'])**(dkt.s))*(Pfac**(dkt.t))*(Afac**(dkt.u))*(z**(dkt.v))    #Ecuacion de regrecion KT
            tabla.loc[i,'KT']=dkt['computkt'].sum()                                                                 
            
            dkq['computkq']=(dkq.Tstuv)*((tabla.loc[i,'J'])**(dkq.s))*(Pfac**(dkq.t))*(Afac**(dkq.u))*(z**(dkq.v)) #Ecuacion de regrecion KQ
            tabla.loc[i,'KQ']=dkq['computkq'].sum()
            
            if (re>1e6):                                                                                              #Correcion en funcion al numero de Reynold
                
                delkt=0.000353485                                                               \
                    -0.00333758*(Afac)*((tabla.loc[i,'J'])**2)                                 \
                    -0.00478125*(Afac)*(Pfac)*((tabla.loc[i,'J']))                             \
                    +0.000257792*((log10(re)-0.301)**2)*(Afac)*((tabla.loc[i,'J'])**2)         \
                    +0.0000643192*(log10(re)-0.301)*(Pfac**6)*((tabla.loc[i,'J'])**2)          \
                    -0.0000110636*((log10(re)-0.301)**2)*(Pfac**6)*((tabla.loc[i,'J'])**2)     \
                    -0.0000276305*((log10(re)-0.301)**2)*(z)*(Afac)*((tabla.loc[i,'J'])**2)    \
                    +0.0000954*(log10(re)-0.301)*(z)*(Afac)*(Pfac)*(tabla.loc[i,'J'])          \
                    +0.0000032049*(log10(re)-0.301)*(z**2)*(Afac)*(Pfac**3)*(tabla.loc[i,'J']) 
                
                tabla.loc[i,'KT']=tabla.loc[i,'KT']+delkt
                
                delkq=-0.000591412                                                               \
                    +0.00696898*(Pfac)                                                         \
                    -0.0000666654*(z)*(Pfac**6)                                                \
                    +0.0160818*(Afac**2)                                                       \
                    -0.000938091*(log10(re)-0.301)*(Pfac)                                      \
                    -0.00059593*(log10(re)-0.301)*(Pfac**2)                                    \
                    +0.0000782099*((log10(re)-0.301)**2)*(Pfac**2)                             \
                    +0.0000052199*(log10(re)-0.301)*(z)*(Afac)*((tabla.loc[i,'J'])**2)         \
                    -0.00000088528*((log10(re)-0.301)**2)*(z)*(Afac)*(Pfac)*(tabla.loc[i,'J']) \
                    +0.0000230171*(log10(re)-0.301)*(z)*(Pfac**6)                              \
                    -0.00000184341*((log10(re)-0.301)**2)*(z)*(Pfac**6)                        \
                    -0.00400252*(log10(re)-0.301)*(Afac**2)                                    \
                    +0.000220915*((log10(re)-0.301)**2)*(Afac**2)
                
                tabla.loc[i,'KQ']=tabla.loc[i,'KQ']+(delkq)
        
        tabla['EFIC']=((tabla.J)*(tabla.KT))/(2*pi*tabla.KQ) #ARREGLAR EFICIENCIA
        return tabla


#Datos Ejemplo
Pfac=1.08
Afac=0.674204932
z=4
re=1e7

KQKT=SerieB(Pfac,Afac,z,re)
tabla=KQKT.Curves_one(Pfac,Afac,z,re)
print(tabla)
tabla.plot(x="J",y=["KT","KQ"])
plt.show()


#KQKT.Write_Datos()






#tabla.plot(x='J',y=['KT','10KQ'])
#tabla.plot(x="J",y=["KT","KQ"])
#tabla.plot(x="J",y=["KT","10KQ"],kind="scatter" ) #kind="scatter"
#tabla.plot("J","10KQ",marker='+',color='darkblue',linestyle='-.',kind="scatter") #kind="scatter"

#plt.show()