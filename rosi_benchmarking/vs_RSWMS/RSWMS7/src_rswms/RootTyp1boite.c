/*****************************************************************************/
/*     ROOT TYP, 09/03/02                                                    */
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "RootTyp1boite.h"

const int DeltaT=2;         /* Pas de temps, en jours */
const double Epsilon=1.0e-6; /* Petite valeur, proche de 0 */
const double Pi=3.14159265;  /* Valeur approchée de la constante Pi */
const double EpaissHor=5.0;  /* Epaisseur horizons de sol */

/* Fichiers */
FILE *FNoeud;  /* Contient la structure, sous forme de noeuds */
FILE *FPar;    /* Paramètres */
FILE *FSol;    /* Informations sur le sol, par horizons */
FILE *FLim;    /* Informations sur les limites du sol

/* Paramètres, lus dans fichier paramètres */
int P_Duree; /* Durée de la simulation */
int P_NbVag; /* Nombre de vagues de réitération */
int P_NbAxeOG; /* Nombre d'axes au départ sur l'organe générateur */
int P_TpsReitVag[NBVAGMAX+1]; /* Vagues de réitération */
float P_CoeffCroissRad;
float P_AngIMoy[NBTYPMAX]; /* Angle d'insertion ramification moyenne */
float P_AngIEt[NBTYPMAX];  /* Angle d'insertion ramification écart-type */
float P_AngIReitMoy[NBTYPMAX]; /* Angle d'insertion réitération moyenne */
float P_AngIReitEt[NBTYPMAX];  /* Angle d'insertion réitération écart-type */
float P_DurDevPrim[NBTYPMAX];  /* Durée de développement du primordium */
float P_CroissMoy[2][NBTYPMAX]; /* Paramètres de la courbe de croissance (moy) */
float P_CroissEt[2][NBTYPMAX]; /* Paramètres de la courbe de croissance (e-t) */
float P_RamifMoy[NBTYPMAX]; /* Distance inter-ramification (moyenne) */
float P_RamifEt[NBTYPMAX]; /* Distance inter-ramification (écart-type)*/
int P_TTrop[NBTYPMAX]; /* Type de tropisme (0: plagio; -1: geo-; +1: geo+; 2: exo */
float P_ITrop[NBTYPMAX]; /* Intensité du tropisme */
float P_SCMeca[NBTYPMAX]; /* Sensibilité contrainte mécanique */
float P_PropTypRamif[NBTYPMAX][NBTYPMAX]; /* Proportion des types de ramif */
float P_DureeNecrose[NBTYPMAX]; /* Durée entre arrêt de croissance et nécrose */
float P_DiamPrim[NBTYPMAX]; /* Diamètre primaire (de la pointe) */
float P_ProbReit[NBTYPMAX]; /* Probabilité de réitération durant le pas (quand vague) */
int P_NbReitMin[NBTYPMAX]; /* Nombre minimal de réitérations, lorsque réitération */
int P_NbReitMax[NBTYPMAX]; /* Nombre maximal de réitérations, lorsque réitération */
float P_AgeTransf[NBTYPMAX]; /* Age de déclenchement possibilité de transformation */
float P_ProbTransf[NBTYPMAX]; /* Probabilité de transformation durant une journée */
int P_SensTransf[NBTYPMAX]; /* Sens de la transformation : -1 ou +1 */
int P_TypeDest[NBTYPMAX]; /* Type à faire apparaître dans le fichier out noeud */

/* Variables limites du système*/
float Xmin;
float Xmax;
float Ymin;
float Ymax;
float Zmax;

/* Variables globales diverses */
int Temps=0;  /* Le temps, en jours */
r3 Orig; /* Position d'origine de l'organe générateur */

PTSysRac SR;  /* Le système racinaire */
TSol Sol;     /* Le sol */

/****************************************************************************/
double FRandUnif(void)
     /* Cette fonction tire un aléatoire uniforme réel entre 0 et 1 */
{
  double tirage;
  tirage=(double) rand()/(double) RAND_MAX;
  if (tirage<Epsilon) { tirage=Epsilon; }
  return(tirage);
}
/****************************************************************************/
/****************************************************************************/
void Norme(r3 u, r3 un)
     /* Cette fonction norme le vecteur u de l'espace de dimension 3.
	Le vecteur norme de retour est appele un. */
{
  double NorU;
  NorU=sqrt((u[0]*u[0])+(u[1]*u[1])+(u[2]*u[2]));
  if (NorU<Epsilon)
    {
      printf("ATTENTION, vecteur nul ! Sa norme vaut : %f \n",NorU);
      exit(1);
    }
  else
    {
      un[0]=u[0]/NorU;
      un[1]=u[1]/NorU;
      un[2]=u[2]/NorU;
    }
}  /* Fonction Norme */
/****************************************************************************/
/****************************************************************************/
double ProdScal(r3 u,r3 v)
     /* Cette fonction retourne le produit scalaire de 2 vecteurs u et v de
	l'espace a 3 dimensions. */
{
  double ProdScal;
  ProdScal=(u[0]*v[0])+(u[1]*v[1])+(u[2]*v[2]);
  return(ProdScal);
}  /* Fonction ProdScal */
/****************************************************************************/
/****************************************************************************/
void ProdVect(r3 u, r3 v, r3 u_vect_v)
     /* Cette fonction calcule le produit vectoriel de deux vecteurs u et v
	de l'espace de dimension 3. Le vecteur resultant est u_vect_v. */
{
  u_vect_v[0]=(u[1]*v[2])-(v[1]*u[2]);
  u_vect_v[1]=(u[2]*v[0])-(v[2]*u[0]);
  u_vect_v[2]=(u[0]*v[1])-(v[0]*u[1]);
}   /* Fonction ProdVect */
/****************************************************************************/
/****************************************************************************/
void RotVect(double omega, r3 u, r3 x, r3 rot_x)

     /* Cette fonction calcule le vecteur rot_x dans l'espace de dimension 3,
	issu de la rotation du vecteur x autour d'un axe dont u est un vecteur
	unitaire. La rotation se fait d'un angle omega radians. Elle appelle
	PRODSCAL, PRODVECT. */
{
  double uscalx;   /* produit scalaire u.x  */
  r3    uvectx;   /* produit vectoriel u^x */

  uscalx=ProdScal(u,x);
  ProdVect(u,x,uvectx);

  rot_x[0]=((1-cos(omega))*uscalx*u[0])
    +(cos(omega)*x[0])+(sin(omega)*uvectx[0]);
  rot_x[1]=((1-cos(omega))*uscalx*u[1])
    +(cos(omega)*x[1])+(sin(omega)*uvectx[1]);
  rot_x[2]=((1-cos(omega))*uscalx*u[2])
    +(cos(omega)*x[2])+(sin(omega)*uvectx[2]);

}  /* Fonction RotVect */
/****************************************************************************/
/****************************************************************************/
void RotZ(r3 u, r3 v, double teta)
     /* Cette fonction fait tourner "u" d'un angle "teta" autour de l'axe (Oz);
	le vecteur calcule est "v" */
{
  v[0]=(u[0]*cos(teta))-(u[1]*sin(teta));
  v[1]=(u[0]*sin(teta))+(u[1]*cos(teta));
  v[2]=u[2];
}
/****************************************************************************/
/****************************************************************************/
int IRandUnif(int imax)

     /* Cette fonction tire un aléatoire uniforme entier entre 0 et imax */
{
  int tirage;

  tirage=imax+1;
  while (tirage>imax) tirage=rand();
  return tirage;
}
/****************************************************************************/
/****************************************************************************/
void OuvreFichiers(void)
     /* Cette fonction ouvre les fichiers */
{
  FNoeud=fopen("in/noeud.txt","w+");
  FPar=fopen("in/param.txt","rt");
  FSol=fopen("in/sol.txt","rt");
  FLim=fopen("in/limites.txt","rt");
} /* Fonction OuvreFichiers */
/****************************************************************************/
/****************************************************************************/
void LitSol(void)
     /* Fonction de lecture des caractéristiques du sol, une ligne par horizon */
{
  int hor;              /* Compteur des horizons */
  char bid[MAXLINE];    /* Chaîne qui accueille les caractères supplémentaires */

  fgets(bid,MAXLINE-1,FSol);          /* Ligne entête */
  for (hor=0; hor<NBHORMAX; hor++)
    {
      fscanf(FSol,"%f %f %f %d",&Sol[hor].Croiss,&Sol[hor].Ramif,&Sol[hor].ICMeca,&Sol[hor].OCMeca);
      /*  fscanf(FSol,"%f",&Sol[hor].Croiss); /* Favorable à la croissance */
      /*  fscanf(FSol,"%f",&Sol[hor].Ramif);  /* Favorable à la ramification */
      /*  fscanf(FSol,"%f",&Sol[hor].ICMeca); /* Intensité de la contrainte */
      /*  fscanf(FSol,"%d",&Sol[hor].OCMeca); /* Orientation 0: iso, 1: verticale */
      fgets(bid,MAXLINE-1,FSol);
    }

} /* Fonction LitSol */
/****************************************************************************/
/****************************************************************************/
void LitLim(void)
     /* Fonction de lecture des limites du lysimètre*/
{
  char bid[MAXLINE];    /* Chaîne qui accueille les caractères supplémentaires */


  fgets(bid,MAXLINE-1,FLim);          /* Ligne entête */
  fscanf(FLim,"%f %f %f %f %f",&Xmin,&Xmax,&Ymin,&Ymax,&Zmax);

}/* Fonction LitLim */
/****************************************************************************/
/****************************************************************************/
double CroissSol(TSol Sol, double Profondeur)
     /* Renvoie le coefficient de croissance du sol à la Profondeur donnée */
{
  int Hor;

  Hor=(int) floor(-Profondeur/EpaissHor);
  if (Hor>=NBHORMAX) Hor=NBHORMAX-1;
  if (Hor<0) Hor=0;

  return(Sol[Hor].Croiss);
} /* Fonction CroissSol */
/****************************************************************************/
/****************************************************************************/
double RamifSol(TSol Sol, double Profondeur)
     /* Renvoie le coefficient de ramification du sol à la Profondeur donnée */
{
  int Hor;

  Hor=(int) floor(-Profondeur/EpaissHor);
  if (Hor>=NBHORMAX) Hor=NBHORMAX-1;
  if (Hor<0) Hor=0;

  return(Sol[Hor].Ramif);
} /* Fonction RamifSol */
/****************************************************************************/
/****************************************************************************/
double ICMecaSol(double Profondeur)
     /* Renvoie l'intensité de la contraine méca du sol à la Profondeur donnée */
{
  int Hor;

  Hor=(int) floor(-Profondeur/EpaissHor);
  if (Hor>=NBHORMAX) Hor=NBHORMAX-1;
  if (Hor<0) Hor=0;

  return(Sol[Hor].ICMeca);
} /* Fonction ICMecaSol */
/****************************************************************************/
/****************************************************************************/
int OCMecaSol(double Profondeur)
     /* Renvoie l'indice de la direction de contrainte : 0 pour iso, 1 pour verti */
{
  int Hor;

  Hor=(int) floor(-Profondeur/EpaissHor);
  if (Hor>=NBHORMAX) Hor=NBHORMAX-1;
  if (Hor<0) Hor=0;

  return(Sol[Hor].OCMeca);
} /* Fonction OCMecaSol */
/****************************************************************************/
/****************************************************************************/
double TireAngI(int TypeFils)
{  /* Tire l'angle d'insertion d'une ramification sur sa mère */
  double TireAngI,tire1,tire2;

  tire1=FRandUnif();
  tire2=FRandUnif();
  TireAngI=P_AngIMoy[TypeFils]+(P_AngIEt[TypeFils]
				*sqrt(-log(tire1))*cos(Pi*tire2)*1.414);
  return(TireAngI);
} /* Fonction TireAngI */
/****************************************************************************/
/****************************************************************************/
double TireAngIReit(int Type)
{    /* Tire l'angle d'insertion d'une réitération sur sa mère */
  double TireAngIReit,tire1,tire2;

  tire1=FRandUnif();
  tire2=FRandUnif();
  TireAngIReit=P_AngIReitMoy[Type]+(P_AngIReitEt[Type]*sqrt(-log(tire1))*cos(Pi*tire2)*1.414);
  return(TireAngIReit);
} /* Fonction TireAngIReit */
/****************************************************************************/
/****************************************************************************/
double TireAngRad(void)
{   /* Tire l'angle radial dans l'intervalle 0 - 2*Pi */

  return (2.0*Pi*FRandUnif());
} /* Fonction TireAngRad */
/****************************************************************************/
/****************************************************************************/
void IncreNbNoeudSR(PTSysRac SR)
     /* Incrémente le nombre de noeuds qui a été formé dans ce système SR */
{
  SR->NbNoeudForm++;
} /* Fonction IncreNbNoeudSR */
/****************************************************************************/
/****************************************************************************/
PTNoeud CreeNoeud(void)
     /* Cette fonction retourne une nouvelle variable de type PTNoeud,
	c'est-à-dire un pointeur sur le type TNoeud */
{
  PTNoeud Nd;
  Nd=(PTNoeud) malloc(sizeof(TNoeud));
  if (Nd==NULL)
    { printf("Problème mémoire allocation noeud \n"); exit(1); }

  return Nd;
} /* Fonction CreeNoeud */
/****************************************************************************/
/****************************************************************************/
PTNoeud InitialiseNoeud(long int Num, r3 Position, double Diam, PTAxe Pere, PTAxe Fils)
     /* Cette fonction retourne une nouvelle variable de type PTNoeud,
	dont une partie des valeurs est initialisée */
{
  PTNoeud Nd;

  Nd=CreeNoeud();

  Nd->Num=Num;
  Nd->JourForm=Temps;
  Nd->Necrose=0;

  Nd->Diametre=Diam;
  Nd->Pere=Pere;
  Nd->Fils=Fils;

  Nd->Pos[0]=Position[0];
  Nd->Pos[1]=Position[1];
  Nd->Pos[2]=Position[2];

  Nd->SuivSPere=NULL;
  Nd->SuivSFils=NULL;
  Nd->Prec=NULL;
  Nd->SuivSAxe=NULL;

  return Nd;
} /* Fonction InitialiseNoeud */
/****************************************************************************/
/****************************************************************************/
void DetruitNoeud(PTNoeud NdADetruire)
     /* Supprime un noeud en mémoire */
{
  free(NdADetruire);
} /* Fonction DetruitNoeud */
/****************************************************************************/
/****************************************************************************/
PTMeristeme CreeMeris(void)
     /* Cette fonction retourne une nouvelle variable de type PTMeristeme,
	c'est-à-dire un pointeur sur le type TMeristeme */
{
  PTMeristeme Meris;
  Meris=(PTMeristeme) malloc(sizeof(TMeristeme));
  if (Meris==NULL)
    { printf("Problème mémoire allocation Meristème \n"); exit(1); }

  return Meris;
} /* Fonction CreeMeris */
/****************************************************************************/
/****************************************************************************/
void TirePCroissMeris(PTMeristeme Meris)
{
  double tire1,tire2;
  const double PCroissMin=0.3;

  tire1=FRandUnif();
  tire2=FRandUnif();

  /* On tire une asymptote suivant une loi normale et une vitesse initiale fixe */
  Meris->PCroiss[0]=P_CroissMoy[0][Meris->Type]+(P_CroissEt[0][Meris->Type]*
						 sqrt(-log(tire1))*cos(Pi*tire2)*1.414);  /* Asymptote */
  if (Meris->PCroiss[0]<PCroissMin) { Meris->PCroiss[0]=PCroissMin; };

  Meris->PCroiss[1]=P_CroissMoy[1][Meris->Type]; /* Vitesse initiale */
} /* Fonction TirePCroissMeris */
/****************************************************************************/
/****************************************************************************/
void TirePRamifMeris(PTMeristeme Meris)
{  /* Affectation de la distance inter-ramification */
  if (Meris->PCroiss[0]>0.01) { Meris->PRamif=P_RamifMoy[Meris->Type]; }
  else   { Meris->PRamif=1000.0; }
} /* Fonction TirePRamifMeristeme */
/****************************************************************************/
/****************************************************************************/
PTMeristeme InitialiseMeris(int Type, r3 Position, r3 Direction)
     /* Cette fonction retourne une nouvelle variable de type PTMeristeme,
	dont les valeurs sont en partie initialisées */
{
  PTMeristeme Meris;

  Meris=CreeMeris();

  Meris->DistPrimInit=0.0;
  Meris->Age=0.0;
  Meris->Type=Type;
  Meris->TypeDest=P_TypeDest[Meris->Type];
  Meris->Diametre=P_DiamPrim[Meris->Type];
  Meris->Senile=0;
  Meris->Mature=0;

  Meris->Coord[0]=Position[0];
  Meris->Coord[1]=Position[1];
  Meris->Coord[2]=Position[2];

  Meris->DirCroiss[0]=Direction[0];
  Meris->DirCroiss[1]=Direction[1];
  Meris->DirCroiss[2]=Direction[2];

  Meris->DirInit[0]=Direction[0];
  Meris->DirInit[1]=Direction[1];
  Meris->DirInit[2]=Direction[2];

  TirePCroissMeris(Meris);
  TirePRamifMeris(Meris);

  return Meris;
} /* Fonction InitialiseMeristeme */
/****************************************************************************/
/****************************************************************************/
void DeflecMecaMeris(PTMeristeme Meris, r3 DirApresMeca, double Elong)
{
  const double Teta=15.0; /* Angle autour de G, en degres */

  r3 VTire,VTireN,DirInt;
  double Profondeur, Cont, Aa, Rr, X, Y; /* Aa et Rr rajoutés suite interv. Xavier */

  X=Meris->Coord[0];
  Y=Meris->Coord[1];
  Profondeur=Meris->Coord[2];
  Cont=ICMecaSol(Profondeur);
  if (OCMecaSol(Profondeur)==1) /* Contrainte anisotrope verticale */
    {
      /* Tirage vecteur dans l'angle Teta autour de G */
      Aa=FRandUnif()*2*Pi; /* Angle aléatoire */
      Rr=sqrt(FRandUnif())*sin(Pi*Teta/180);
      VTireN[0]=Rr*cos(Aa);
      VTireN[1]=Rr*sin(Aa);
      VTireN[2]=-sqrt(1.0-(VTireN[0]*VTireN[0])-(VTireN[1]*VTireN[1]));

      DirInt[0]=Meris->DirCroiss[0]+(Elong*VTireN[0]*Cont*P_SCMeca[Meris->Type]);
      DirInt[1]=Meris->DirCroiss[1]+(Elong*VTireN[1]*Cont*P_SCMeca[Meris->Type]);
      DirInt[2]=Meris->DirCroiss[2]+(Elong*VTireN[2]*Cont*P_SCMeca[Meris->Type]);
    }
  else    /* Contrainte isotrope [OCMecaSol(Profondeur)==0] */
    {
      VTire[0]=2.0*FRandUnif()-1.0;
      VTire[1]=2.0*FRandUnif()-1.0;
      VTire[2]=-2.0*FRandUnif()+1.0;
      Norme(VTire,VTireN);
      if (ProdScal(VTireN,Meris->DirCroiss)<0.0)
	{
	  VTireN[0]=-VTireN[0];
	  VTireN[1]=-VTireN[1];
	  VTireN[2]=-VTireN[2];
	}
      DirInt[0]=Meris->DirCroiss[0]+(Elong*VTireN[0]*Cont*P_SCMeca[Meris->Type]);
      DirInt[1]=Meris->DirCroiss[1]+(Elong*VTireN[1]*Cont*P_SCMeca[Meris->Type]);
      DirInt[2]=Meris->DirCroiss[2]+(Elong*VTireN[2]*Cont*P_SCMeca[Meris->Type]);
    }

  /*javaux new*/
  if (abs(Profondeur)>=abs(Zmax))
     {
      Cont=100000;
     
      Aa=FRandUnif()*Pi; /* Angle aléatoire */
      Rr=1-sin(Pi*5/180);

      if (X>=Xmax)
        {
        Aa=Aa+Pi/2;
        }
      else if(X<=Xmin)
        {
        Aa=Aa-Pi/2;
        }

      if ((Y>=Ymax))
        {
         Aa=Aa+Pi;
        }
 
      VTireN[0]=Rr*cos(Aa);
      VTireN[1]=Rr*sin(Aa);
      VTireN[2]=sqrt(1.0-(VTireN[0]*VTireN[0])-(VTireN[1]*VTireN[1]));

      DirInt[0]=Meris->DirCroiss[0]+(Elong*VTireN[0]*Cont);
      DirInt[1]=Meris->DirCroiss[1]+(Elong*VTireN[1]*Cont);
      DirInt[2]=Meris->DirCroiss[2]+(Elong*VTireN[2]*Cont);

      }

  else if ((X>=Xmax)|(X<=Xmin)|(Y>=Ymax)|(Y<=Ymin)) /*hors limite*/
     {
/*     printf("X,Y,Z :");
        printf("% 10.3f % 10.3f % 10.3f",X,Y,Profondeur);
        printf("minX,maxX,minY,maxY,maxZ :");
        printf("% 10.3f % 10.3f % 10.3f % 10.3f % 10.3f",Xmin,Xmax,Ymin,Ymax,Zmax);
       Tirage vecteur dans l'angle Teta autour de G */
      Cont=100000;
      Aa=sqrt(FRandUnif())*Pi; /* Angle aléatoire */


      if (X>=Xmax)
        {
        Aa=Aa+Pi/2;
        }
      else if(X<=Xmin)
        {
        Aa=Aa-Pi/2;
        }

      if ((Y>=Ymax))
        {
         Aa=Aa+Pi;
        }
      else if(Y<=Ymin)
        {
        
        }  

      Rr=(FRandUnif())*sin(Pi*25/180);
      VTireN[0]=Rr*cos(Aa);
      VTireN[1]=Rr*sin(Aa);
      VTireN[2]=-sqrt(1.0-(VTireN[0]*VTireN[0])-(VTireN[1]*VTireN[1]));

      DirInt[0]=Meris->DirCroiss[0]+(Elong*VTireN[0]*Cont);
      DirInt[1]=Meris->DirCroiss[1]+(Elong*VTireN[1]*Cont);
      DirInt[2]=Meris->DirCroiss[2]+(Elong*VTireN[2]*Cont);
     }

  Norme(DirInt,DirApresMeca);


} /* Fonction DeflecMecaMeris */
/****************************************************************************/
/****************************************************************************/
void DeflecGeoMeris(PTMeristeme Meris, r3 DirApresMeca, r3 DirApresGeo, double Elong)
     /* Version avec plagiotropisme */
{
  r3 DirInt,VGeoInt,VGeo;

  switch (P_TTrop[Meris->Type]) {
  case -1 : VGeo[0]=0.0;                  /* Gravitropisme négatif */
    VGeo[1]=0.0;
    VGeo[2]=1.0;
    break;
  case 0 : VGeoInt[0]=Meris->DirInit[0]; /* Plagiotropisme */
    VGeoInt[1]=Meris->DirInit[1];
    VGeoInt[2]=0.0;
    Norme(VGeoInt,VGeo);
    break;
  case 1 : VGeo[0]=0.0;                  /* Gravitropisme positif */
    VGeo[1]=0.0;
    VGeo[2]=-1.0;
    break;
  case 2 : VGeoInt[0]=Meris->DirInit[0]; /* Exotropisme */
    VGeoInt[1]=Meris->DirInit[1];
    VGeoInt[2]=Meris->DirInit[2];
    Norme(VGeoInt,VGeo);
    break;
  default : VGeo[0]=0.0;                 /* Gravitropisme positif */
    VGeo[1]=0.0;
    VGeo[2]=-1.0;
    break;
  }

  DirInt[0]=DirApresMeca[0]+(VGeo[0]*P_ITrop[Meris->Type]*Elong);
  DirInt[1]=DirApresMeca[1]+(VGeo[1]*P_ITrop[Meris->Type]*Elong);
  DirInt[2]=DirApresMeca[2]+(VGeo[2]*P_ITrop[Meris->Type]*Elong);

  Norme(DirInt,DirApresGeo);
} /* Fonction DeflecGeoMeris */
/****************************************************************************/
/****************************************************************************/
void DeflecSurfMeris(PTMeristeme Meris, r3 DirApresGeo, r3 DirApresSurf)
{
  const double ProfLim=-3.0*FRandUnif();
  r3 DirInt;
  DirInt[0]=DirApresGeo[0];
  DirInt[1]=DirApresGeo[1];
  DirInt[2]=DirApresGeo[2];

  if ((DirInt[2]>0.0) && ((Meris->Coord[2])>ProfLim) && (Meris->Type>0))
    DirInt[2]=DirInt[2]/8.0;
  Norme(DirInt,DirApresSurf);
} /* Fonction DeflecSurfMeris */
/****************************************************************************/
/****************************************************************************/
void ReorienteMeris(PTMeristeme Meris, double Elong)
{
  r3 DirInt1, DirInt2, NouvDir;

  DeflecMecaMeris(Meris,DirInt1,Elong);
  DeflecGeoMeris(Meris,DirInt1,DirInt2,Elong);
  DeflecSurfMeris(Meris,DirInt2,NouvDir);

  Meris->DirCroiss[0]=NouvDir[0];
  Meris->DirCroiss[1]=NouvDir[1];
  Meris->DirCroiss[2]=NouvDir[2];


} /* Fonction ReorienteMeris */
/****************************************************************************/
/****************************************************************************/
double CalcElongationMeris(PTMeristeme Meris)

{
  /* Calcul de l'élongation potentielle en subdivisant l'intervalle de temps,
     et prise en compte du coefficient de croissance du sol à la profondeur ad hoc */
  int NbDiv=20;
  double dt,t,A,b,ElongPot;
  /* int i; */

  dt=(double) DeltaT/(double) NbDiv;
  t=Meris->Age; /* Age méristème vrai */
  A=Meris->PCroiss[0];  /* Asymptote */
  b=Meris->PCroiss[1];  /* Vitesse initiale */
  ElongPot = -A * (exp(-b * (t + DeltaT) / A) - exp(-b * t / A));
  /*  ElongPot=0.0;
      for (i=1; i<=NbDiv; i++)
      { ElongPot=ElongPot+(b*exp(-b*(t+(i*dt))/A)*dt); } */
  return(ElongPot*CroissSol(Sol,Meris->Coord[2]));

} /* Fonction CalcElongationMeris */
/****************************************************************************/
/****************************************************************************/
void VieillitMeris(PTMeristeme Meris)
{ /* Incrémente l'âge du méristème selon le pas de temps */
  Meris->Age=Meris->Age+DeltaT;
} /* Fonction VieillitMeris */
/****************************************************************************/
/****************************************************************************/
void MatureMeris(PTMeristeme Meris)
{ /* Assure l'évolution du primordium en méristème si son âge est atteint */
  if ((!Meris->Mature)&&(Meris->Age>P_DurDevPrim[Meris->Type]))
    {
      Meris->Mature=1;  /* Le primordium devient méristème vrai */
      Meris->Age=0.0;   /* Son âge est réinitialisé à 0 en tant que méristème */
    }
} /* Fonction MatureMeris */
/****************************************************************************/
/****************************************************************************/
void SenesceMeris(PTMeristeme Meris)
{ /* Rend sénescent le méristème qui ne s'allonge plus ou qui a réitéré */
  Meris->Senile=1;
} /* Fonction SenesceMeris */
/****************************************************************************/
/****************************************************************************/
void TransformeMeris(PTMeristeme Meris)
{ /* Réalise éventuellement la transformation du type du méristème */
  double ProbaDeTransformation;


  /* Calcule la probabilité de transformation */
  ProbaDeTransformation=((double) Meris->Mature)*
    ((double) Meris->Age>P_AgeTransf[Meris->Type])
    *P_ProbTransf[Meris->Type];
  /*  printf("Proba de transformation : %7.3d", ProbaDeTransformation, "\n"); */

  if (FRandUnif()<ProbaDeTransformation) /* transformation */
    {
      Meris->Type=Meris->Type+P_SensTransf[Meris->Type];
      Meris->TypeDest=P_TypeDest[Meris->Type];
      Meris->Age=0.0;  /* Age réinitialisé suite à transformation */
      TirePCroissMeris(Meris);
      TirePRamifMeris(Meris);
      Meris->DirInit[0]=Meris->DirCroiss[0];
      Meris->DirInit[1]=Meris->DirCroiss[1];
      Meris->DirInit[2]=Meris->DirCroiss[2];
    }
  if ((Meris->Type<0)||(Meris->Type>NBTYPMAX))
    {
      printf("Problème dans TransformeMeris, Type du méristème non conforme \n");
      exit(1);
    }
} /* Fonction TransformeMeris */
/****************************************************************************/
/****************************************************************************/
void DeveloppeMeris(PTMeristeme Meris)
{ /* Assure l'évolution du méristème */
  VieillitMeris(Meris);
  MatureMeris(Meris);
  TransformeMeris(Meris);
} /* Fonction DeveloppeMeris */
/****************************************************************************/
/****************************************************************************/
void DeplaceMeris(PTMeristeme Meris, double Elong)
{ /* Assure le déplacement du méristème suite à croissance axiale */

  /* Sa position est modifiée */
  Meris->Coord[0]=Meris->Coord[0]+(Elong*Meris->DirCroiss[0]);
  Meris->Coord[1]=Meris->Coord[1]+(Elong*Meris->DirCroiss[1]);
  Meris->Coord[2]=Meris->Coord[2]+(Elong*Meris->DirCroiss[2]);

  /* Son attribut DistPrimInit est modifié */
  Meris->DistPrimInit=Meris->DistPrimInit+Elong;

} /* Fonction DeplaceMeris */
/****************************************************************************/
/****************************************************************************/
double DistInterRamifMeris(PTMeristeme Meris, TSol Sol)
{ /* Renvoie la valeur locale de la distance inter-ramification du méristème */

  return (Meris->PRamif*RamifSol(Sol,Meris->Coord[2]));

} /* Fonction DistInterRamifMeris */
/****************************************************************************/
/****************************************************************************/
void DetruitMeris(PTMeristeme MerisADetruire)
     /* Supprime un méristème */
{
  free(MerisADetruire);
} /* Fonction DetruitMeris
     /****************************************************************************/
/****************************************************************************/
PTAxe CreeAxe(void)
     /* Cette fonction retourne une nouvelle variable de type PTAxe,
	c'est-à-dire un pointeur sur le type TAxe */
{
  PTAxe Axe;
  Axe=(PTAxe) malloc(sizeof(TAxe));
  if (Axe==NULL)
    { printf("Problème mémoire allocation dans CreeAxe \n"); exit(1); }

  return Axe;
} /* Fonction CreeAxe */
/****************************************************************************/
/****************************************************************************/
PTAxe InitialiseAxe(long int NumAxe, int TypeMeris, r3 Origine, r3 DirInit, PTAxe AxePere)
     /* Cette fonction retourne une nouvelle variable de type PTAxe,
	c'est-à-dire un pointeur sur le type TAxe */
{
  PTAxe NouvAxe;
  PTMeristeme Meris;
  PTNoeud PremierNoeud;

  NouvAxe=CreeAxe();
  Meris=InitialiseMeris(TypeMeris,Origine,DirInit);
  PremierNoeud=InitialiseNoeud(SR->NbNoeudForm+1,Origine,P_DiamPrim[TypeMeris],AxePere,NouvAxe);
  NouvAxe->Meris=Meris;
  NouvAxe->PremNoeud=PremierNoeud;
  NouvAxe->DernNoeud=PremierNoeud;
  NouvAxe->NbNoeud=1;
  NouvAxe->Num=NumAxe;
  NouvAxe->Pere=AxePere;

  NouvAxe->Suivant=NULL;
  NouvAxe->Precedent=NULL;

  return NouvAxe;
} /* Fonction InitialiseAxe */
/****************************************************************************/
/****************************************************************************/
PTNoeud AvtDernNoeudAxe(PTAxe Axe)
     /* Cette fonction retourne une variable de type PTNoeud,
	qui pointe sur l'avant dernier noeud de l'axe, s'il existe */
{
  PTNoeud NdCour,NdPrec;
  NdCour=Axe->PremNoeud;
  if (NdCour->SuivSFils==NULL)
    { printf("Pas d'avant-dernier noeud, car un seul noeud \n"); exit(1); }
  else
    {
      NdPrec=NdCour;
      NdCour=NdCour->SuivSFils;
      while ((NdCour->SuivSFils!=NULL)||(NdCour->SuivSPere!=NULL))
	{
	  NdPrec=NdCour;
	  if (NdCour->SuivSPere!=NULL) NdCour=NdCour->SuivSPere;
	  else NdCour=NdCour->SuivSFils;
	}
    }
  return NdPrec;
} /* Fonction AvtDernNoeudAxe */
/****************************************************************************/
/****************************************************************************/
void AjouteNoeudTermAxe(PTAxe Axe, PTNoeud NdAAjouter)
     /* Cette fonction ajoute un noeud en position terminale (apicale)
	à l'axe concerné, et incrémente son compteur de noeuds */
{
  PTNoeud AncienNdTerm;

  AncienNdTerm=Axe->DernNoeud;

  AncienNdTerm->SuivSFils=NdAAjouter;
  if (AncienNdTerm->SuivSAxe==NULL)
      AncienNdTerm->SuivSAxe=NdAAjouter;
  NdAAjouter->Prec=AncienNdTerm;
  NdAAjouter->Pere=NULL;
  Axe->DernNoeud=NdAAjouter;
  Axe->NbNoeud++;

} /* Fonction AjouteNoeudTermAxe */
/****************************************************************************/
/****************************************************************************/
void AjouteNoeudLatAxe(PTAxe Axe, PTNoeud NdAAjouter)
     /* Cette fonction ajoute un noeud en position latérale (lors d'une ramification
	ou d'une réitération) de l'axe concerné, et incrémente son compteur de noeuds */
{
  PTNoeud NdPrec;

  NdPrec=Axe->DernNoeud->Prec;
  if (!(NdPrec->SuivSPere!=NULL && NdPrec->SuivSFils!=NULL))
      NdPrec->SuivSAxe=NdAAjouter;
  NdAAjouter->Prec=NdPrec;
  NdAAjouter->SuivSPere=Axe->DernNoeud;
  NdAAjouter->SuivSFils=NULL;
  NdAAjouter->SuivSAxe=Axe->DernNoeud;

  NdAAjouter->Pere=Axe;
  Axe->NbNoeud++;

  /* Reaffecter la succession du precedent */
  /* Si le précédent noeud était un noeud de branchement */
  if (NdPrec->SuivSPere == NdAAjouter->SuivSPere) NdPrec->SuivSPere=NdAAjouter;

  /* Si le précédent noeud n'était pas un noeud de branchement */
  if (NdPrec->SuivSFils == NdAAjouter->SuivSPere) NdPrec->SuivSFils=NdAAjouter;

  Axe->DernNoeud->Prec=NdAAjouter;

} /* Fonction AjouteNoeudLatAxe */
/****************************************************************************/
/****************************************************************************/
PTNoeud DernNoeudAxe(PTAxe Axe)
     /* Cette fonction retourne une variable de type PTNoeud,
	qui pointe sur le dernier noeud de l'axe, le plus distal */
{
  PTNoeud NdCour;
  NdCour=Axe->PremNoeud;
  if (NdCour->SuivSFils!=NULL)
    {
      NdCour=NdCour->SuivSFils;
      while ((NdCour->SuivSFils!=NULL)||(NdCour->SuivSPere!=NULL))
	{
	  if (NdCour->SuivSPere!=NULL) NdCour=NdCour->SuivSPere;
	  else NdCour=NdCour->SuivSFils;
	}
    }
  return NdCour;
} /* Fonction DernNoeudAxe */
/****************************************************************************/
/****************************************************************************/
void DeveloppeAxe(PTAxe Axe)
{
  DeveloppeMeris(Axe->Meris);
} /* Fonction DeveloppeAxe */
/****************************************************************************/
/****************************************************************************/
void AllongeAxe(PTAxe Axe)
{
  const double LongSeuilCroiss=1.0e-2;
  double Elongation;
  PTNoeud NouvNd;

  if ((!Axe->Meris->Senile) && (Axe->Meris->Mature))
    {
      Elongation=CalcElongationMeris(Axe->Meris);
      if (Elongation<LongSeuilCroiss) { SenesceMeris(Axe->Meris); }
      else
	{
	  /* Calcule et affecte la nouvelle direction de croissance du méristème */
	  ReorienteMeris(Axe->Meris,Elongation);

	  /* Le méristème se déplace */
	  DeplaceMeris(Axe->Meris,Elongation);

	  /* Il génère un nouveau noeud sur cet axe à sa nouvelle position */
	  IncreNbNoeudSR(SR);

	  NouvNd=InitialiseNoeud(SR->NbNoeudForm,Axe->Meris->Coord,Axe->Meris->Diametre,NULL,Axe);

	  AjouteNoeudTermAxe(Axe,NouvNd);

	}
    }
} /* Fonction AllongeAxe */
/****************************************************************************/
/****************************************************************************/
void DetruitAxe(PTAxe AxeADetruire)
     /* Supprime un axe en supprimant ses noeuds, puis l'axe lui-même */
{
  PTNoeud NdCour, NdAEnlever;

  /* Liberer tous les noeuds de cet axe */
  NdCour=AxeADetruire->PremNoeud;
  while (NdCour->SuivSFils!=NULL)
    {
      NdAEnlever=NdCour;
      NdCour=NdCour->SuivSFils;
      if (NdCour->SuivSPere!=NULL) { printf("Problème : Axe ramifie a enlever\n"); exit(1); }
      DetruitNoeud(NdAEnlever);
    }
  DetruitNoeud(NdCour); /* Enleve le noeud apical */

  DetruitMeris(AxeADetruire->Meris);

  /* Enlever l'axe en mémoire */
  free(AxeADetruire);

} /* Fonction DetruitAxe */
/****************************************************************************/
/****************************************************************************/
int AxeToutNecrose(PTAxe Axe)
     /* Cette fonction retourne la valeur 1 si l'axe
	a tous ses noeuds necroses et 0 sinon */
{
  PTNoeud NdCour;
  int Resu;

  NdCour=Axe->PremNoeud;
  if (NdCour->SuivSFils==NULL) { Resu=NdCour->Necrose; }   /* Un seul noeud */
  else
    {
      NdCour=NdCour->SuivSFils;
      Resu=NdCour->Necrose;
    }
  return Resu;
} /* Fonction AxeToutNecrose */
/****************************************************************************/
/****************************************************************************/
int AxeToutNecroseAncien(PTAxe Axe)
     /* Cette fonction retourne la valeur 1 si l'axe
	a tous ses noeuds necroses et 0 sinon */
{
  PTNoeud NdCour;
  int Resu=1;   /* Initialisation a la valeur OUI */

  NdCour=Axe->PremNoeud;
  if (NdCour->SuivSFils==NULL) { Resu=NdCour->Necrose; }
  else
    {
      NdCour=NdCour->SuivSFils;
      if (NdCour->Necrose==0) Resu=0;
      while ((NdCour->SuivSFils!=NULL)||(NdCour->SuivSPere!=NULL))
	{
	  if (NdCour->SuivSPere!=NULL) NdCour=NdCour->SuivSPere;
	  else NdCour=NdCour->SuivSFils;
	  if (NdCour->Necrose==0) Resu=0;
	}
    }
  return Resu;
} /* Fonction AxeToutNecroseAncien */
/****************************************************************************/
/****************************************************************************/
void AffecValNecroseAxe(PTAxe Axe, int ValNecrose)
     /* Cette fonction affecte a chacun des noeuds de l'axe
	la valeur de necrose (0 ou 1) */
{
  PTNoeud NdCour;

  NdCour=Axe->PremNoeud;
  NdCour->Necrose=ValNecrose;
  if (NdCour->SuivSFils!=NULL)
    {
      NdCour=NdCour->SuivSFils;
      NdCour->Necrose=ValNecrose;
      while ((NdCour->SuivSFils!=NULL)||(NdCour->SuivSPere!=NULL))
	{
	  if (NdCour->SuivSPere!=NULL) NdCour=NdCour->SuivSPere;
	  else NdCour=NdCour->SuivSFils;
	  NdCour->Necrose=ValNecrose;
	}
    }
} /* Fonction AffecValNecroseAxe */
/****************************************************************************/
/****************************************************************************/
void AffecValNecroseAmont(PTAxe Axe, int ValNecrose)
     /* Cette fonction affecte a chacun des noeuds en amont de l'axe
	la valeur de necrose (0 ou 1) */
{
  PTNoeud NdCour;

  NdCour=Axe->PremNoeud->Prec;
  while (NdCour!=NULL)
    {
      NdCour->Necrose=ValNecrose;
      NdCour=NdCour->Prec;
    }

} /* Fonction AffecValNecroseAmont */
/****************************************************************************/
/****************************************************************************/
void AffecValDiamAxe(PTAxe Axe, float Diam)
     /* Cette fonction affecte a chacun des noeuds de l'axe
	la valeur de diametre Diam */
{
  PTNoeud NdCour;
  NdCour=Axe->PremNoeud;
  NdCour->Diametre=Diam;
  if (NdCour->SuivSFils!=NULL)
    {
      NdCour=NdCour->SuivSFils;
      NdCour->Diametre=Diam;
      while ((NdCour->SuivSFils!=NULL)||(NdCour->SuivSPere!=NULL))
	{
	  if (NdCour->SuivSPere!=NULL) NdCour=NdCour->SuivSPere;
	  else NdCour=NdCour->SuivSFils;
	  NdCour->Diametre=Diam;
	}
    }
} /* Fonction AffecValDiamAxe */
/****************************************************************************/
/****************************************************************************/
void IncreValDiamAmont(PTAxe Axe, double Diam, double Coeff)
     /* Cette fonction incremente le diametre de chacun des noeuds en amont
	de l'axe, en incluant son premier noeud */
{
  PTNoeud NdCour;
  double Section,DiamInit;

  NdCour=Axe->PremNoeud;
  while (NdCour!=NULL)
    {
      DiamInit=NdCour->Diametre;
      Section=(Pi*DiamInit*DiamInit/4.0)+(Pi*Coeff*Diam*Diam/4.0);
      NdCour->Diametre=sqrt(4.0*Section/Pi);
      NdCour=NdCour->Prec;
    }

} /* Fonction IncreValDiamAmont */
/****************************************************************************/
/****************************************************************************/
PTSysRac CreeSR(void)
     /* Cette fonction retourne une nouvelle variable de type PTSysRac,
	c'est-à-dire un pointeur sur le type TSysRac */
{
  PTSysRac SR;
  SR=(PTSysRac) malloc(sizeof(TSysRac));
  if (SR==NULL)
    { printf("Problème mémoire allocation dans CreeSR \n"); exit(1); }

  return SR;
} /* Fonction CreeSR */
/****************************************************************************/
/****************************************************************************/
void AjouteAxeSR(PTSysRac SR, PTAxe AxeAAjouter)
     /* Cette fonction insère un axe dans la chaîne des axes du système racinaire,
	elle incrémente en même temps le compteur d'axes et de noeuds */
{
  if ((SR->NbAxeForm - SR->NbAxeSup)==0)  /* Le système racinaire est vide */
    {
      AxeAAjouter->Suivant=NULL;
      AxeAAjouter->Precedent=NULL;
      SR->PremAxe=AxeAAjouter;
      SR->DernAxe=AxeAAjouter;
    }
  else /* Le système contient déjà des axes, chaînage double */
    {
      AxeAAjouter->Suivant=NULL;
      AxeAAjouter->Precedent=SR->DernAxe;
      SR->DernAxe->Suivant=AxeAAjouter;
      SR->DernAxe=AxeAAjouter;
    }
  SR->NbAxeForm++;
  IncreNbNoeudSR(SR);

} /* Fonction AjouteAxeSR */
/****************************************************************************/
/****************************************************************************/
void EnleveAxeSR(PTSysRac SR, PTAxe AxeAEnlever)
     /* Cette fonction enlève un axe dans la chaîne des axes et libère la memoire */
{
  PTNoeud PrecSPere;

  if ((SR->NbAxeForm - SR->NbAxeSup)==0)  /* Le système racinaire est vide */
    {
      printf("ATTENTION, probleme dans EnleveAxeSR, SR vide \n");
      exit(1);
    }
  else
    {
      SR->NbAxeSup++;

      /* Refaire les chainages dans la liste */
      AxeAEnlever->Precedent->Suivant=AxeAEnlever->Suivant;
      AxeAEnlever->Suivant->Precedent=AxeAEnlever->Precedent;

      /* Refaire les connexions, brancher le noeud precedent sur Pere, au noeud suivant sur Pere */
      if (AxeAEnlever->Pere!=NULL)
	{
	  PrecSPere=AxeAEnlever->PremNoeud->Prec;

	  if (PrecSPere->SuivSPere!=NULL) /* le precedent est un noeud de branchement */
	    {
	      if (PrecSPere->Pere==AxeAEnlever->Pere) { /* le precedent est une ramif du meme axe */
		PrecSPere->SuivSPere=AxeAEnlever->PremNoeud->SuivSPere;
              PrecSPere->SuivSAxe=AxeAEnlever->PremNoeud->SuivSAxe;
              }
	      else /* le precedent est branche sur un autre axe */
		PrecSPere->SuivSFils=AxeAEnlever->PremNoeud->SuivSPere; /* SuivSAxe reste identique */
	    }
	  else {
           PrecSPere->SuivSFils=AxeAEnlever->PremNoeud->SuivSPere;  /* precedent n'est pas noeud branchement */
           PrecSPere->SuivSAxe=AxeAEnlever->PremNoeud->SuivSAxe;
           }

	  AxeAEnlever->PremNoeud->SuivSPere->Prec=PrecSPere;

	}

      DetruitAxe(AxeAEnlever);  /* Detruit ses noeuds, son méristème, et lui-même */
    }
} /* Fonction EnleveAxeSR */
/****************************************************************************/
/****************************************************************************/
PTSysRac InitialiseSR(r3 Origine)
{
  /* Initialisation du système racinaire */

  PTSysRac SR;

  SR=CreeSR();  /* Création d'un système racinaire */

  SR->NbAxeForm=0;  /* Initialisation des variables */
  SR->NbAxeSup=0;
  SR->NbNoeudForm=0;
  SR->PremAxe=NULL;
  SR->DernAxe=NULL;
  SR->NumVReitCourante=0;
  SR->ReitPossible=0;

  SR->Origine[0]=Origine[0];  /* Origine du système racinaire */
  SR->Origine[1]=Origine[1];
  SR->Origine[2]=Origine[2];

  SR->AngDep=2.0*Pi*FRandUnif();  /* Orientation */


  return(SR);
}  /* Fonction InitialiseSR */
/****************************************************************************/
/****************************************************************************/
void InstalleSR(PTSysRac SR)
{
  /* Installation du système racinaire, c'est à dire émission des premiers axes
     de type 0 */

  PTAxe NouvAxe;
  int NumAxeOG;
  r3 VInit, DirInit;
  double AngRot,AngI;

  for ((NumAxeOG=1); (NumAxeOG<=P_NbAxeOG); (NumAxeOG++)) /* Pour tous axes de l'OG */
    {
      /* Calcul de la direction initiale de l'axe */
      AngI=TireAngI(0);
      VInit[0]=sin(AngI);
      VInit[1]=0.0;
      VInit[2]=cos(AngI);
      AngRot=SR->AngDep+(2*Pi*NumAxeOG/P_NbAxeOG);
      RotZ(VInit,DirInit,AngRot);

      /* Génération de l'axe et intégration dans le système racinaire */
      NouvAxe=InitialiseAxe(SR->NbAxeForm+1,0,SR->Origine,DirInit,NULL);
      AjouteAxeSR(SR,NouvAxe);
    }

}  /* Fonction InstalleSR */
/****************************************************************************/
/****************************************************************************/
void EtatReiterationSR(PTSysRac SR)
{
  /* Definit l'état du système racinaire en terme de réitération */

  if (Temps>=P_TpsReitVag[SR->NumVReitCourante+1])
    {
      SR->NumVReitCourante++;
      SR->ReitPossible=1;
    }
  else
    {
      SR->ReitPossible=0;
    }
}  /* Fonction EtatReiterationSR */
/****************************************************************************/
/****************************************************************************/
void LitParam(void)

     /* Fonction de lecture des parametres de la simulation */
{
  int i,typ,t1,t2;
  char bid[MAXLINE];

  for (i=1; i<=4; i++) { fgets(bid,MAXLINE-1,FPar); }
  fscanf(FPar,"%d",&P_Duree);
  fgets(bid,MAXLINE-1,FPar);
  fscanf(FPar,"%d",&P_NbAxeOG);
  fgets(bid,MAXLINE-1,FPar);
  fscanf(FPar,"%d",&P_NbVag);
  fgets(bid,MAXLINE-1,FPar);
  for (i=1; (i<=P_NbVag); i++) { fscanf(FPar,"%d",&P_TpsReitVag[i]); }
  fgets(bid,MAXLINE-1,FPar);
  fscanf(FPar,"%f",&P_CoeffCroissRad);
  fgets(bid,MAXLINE-1,FPar);
  for (typ=0; typ<NBTYPMAX; typ++)
    {
      fgets(bid,MAXLINE-1,FPar);
      t2 = -999; 
      fscanf(FPar,"TYPE %d %d", &t1, &t2);
      P_TypeDest[typ] = (t2 == -999) ? t1 : t2;
      fgets(bid,MAXLINE-1,FPar);
      fscanf(FPar,"%f %f",&P_AngIMoy[typ],&P_AngIEt[typ]);
      fgets(bid,MAXLINE-1,FPar);
      fscanf(FPar,"%f %f",&P_AngIReitMoy[typ],&P_AngIReitEt[typ]);
      fgets(bid,MAXLINE-1,FPar);
      fscanf(FPar,"%f",&P_DurDevPrim[typ]);
      fgets(bid,MAXLINE-1,FPar);
      fscanf(FPar,"%f %f",&P_CroissMoy[0][typ],&P_CroissMoy[1][typ]);
      fgets(bid,MAXLINE-1,FPar);
      fscanf(FPar,"%f %f",&P_CroissEt[0][typ],&P_CroissEt[1][typ]);
      fgets(bid,MAXLINE-1,FPar);
      fscanf(FPar,"%f %f",&P_RamifMoy[typ],&P_RamifEt[typ]);
      fgets(bid,MAXLINE-1,FPar);
      fscanf(FPar,"%d",&P_TTrop[typ]);
      fgets(bid,MAXLINE-1,FPar);
      fscanf(FPar,"%f",&P_ITrop[typ]);
      fgets(bid,MAXLINE-1,FPar);
      fscanf(FPar,"%f",&P_SCMeca[typ]);
      fgets(bid,MAXLINE-1,FPar);
      fscanf(FPar,"%f",&P_DiamPrim[typ]);
      fgets(bid,MAXLINE-1,FPar);
      fscanf(FPar,"%f",&P_DureeNecrose[typ]);
      fgets(bid,MAXLINE-1,FPar);
      fscanf(FPar,"%f",&P_ProbReit[typ]);
      fgets(bid,MAXLINE-1,FPar);
      fscanf(FPar,"%d %d",&P_NbReitMin[typ],&P_NbReitMax[typ]);
      fgets(bid,MAXLINE-1,FPar);
      fscanf(FPar,"%f",&P_AgeTransf[typ]);
      fgets(bid,MAXLINE-1,FPar);
      fscanf(FPar,"%f",&P_ProbTransf[typ]);
      fgets(bid,MAXLINE-1,FPar);
      fscanf(FPar,"%d",&P_SensTransf[typ]);
      fgets(bid,MAXLINE-1,FPar);
      for (i=0; (i<NBTYPMAX); i++) { fscanf(FPar,"%f",&P_PropTypRamif[typ][i]); }
      fgets(bid,MAXLINE-1,FPar);
    }

} /* Fonction LitParam */
/****************************************************************************/
/****************************************************************************/
int AxeRamifiable(PTAxe Axe)
{   /* Renvoie 1 ou 0 suivant que l'axe est ramifiable ou non */

  return(Axe->Meris->DistPrimInit>DistInterRamifMeris(Axe->Meris, Sol));

} /* Fonction AxeRamifiable */
/****************************************************************************/
/****************************************************************************/
int TireTypeFils(PTAxe AxePere)
{   /* Tire le type de la ramification suivant proportions */
  int CompteType;
  double tirage,PropCum;

  tirage=FRandUnif();
  PropCum=0.0;
  CompteType=0;

  while ((tirage>PropCum) && (CompteType<NBTYPMAX))
    {
      CompteType++;
      PropCum=PropCum+P_PropTypRamif[AxePere->Meris->Type][CompteType-1];
    }
  if ((CompteType-1<0)||(CompteType>NBTYPMAX))
    {
      printf("Probleme dans TireTypeFils\n");
      exit(1);
    }
  return(CompteType-1);

} /* Fonction TireTypeFils */
/****************************************************************************/
/****************************************************************************/
void OrigineRamif(PTAxe AxePere, r3 OrigineFils)
{   /* Calcule la position du point d'origine d'une ramification */
  OrigineFils[0]=AxePere->Meris->Coord[0]-
    (AxePere->Meris->DistPrimInit*AxePere->Meris->DirCroiss[0]);
  OrigineFils[1]=AxePere->Meris->Coord[1]-
    (AxePere->Meris->DistPrimInit*AxePere->Meris->DirCroiss[1]);
  OrigineFils[2]=AxePere->Meris->Coord[2]-
    (AxePere->Meris->DistPrimInit*AxePere->Meris->DirCroiss[2]);
} /* Fonction OrigineRamif */
/****************************************************************************/
/****************************************************************************/
void OrienteRamif(PTAxe AxePere, int TypeFils, r3 DirFils)
{   /* Calcule la direction d'une Axe Fils */
  r3 VAxeRot,RotDirCroiss;
  double NorVProjHor,AngRot;

  /* Calcul de la norme de la projection direction sur plan horizontal */
  NorVProjHor=sqrt((AxePere->Meris->DirCroiss[0]*AxePere->Meris->DirCroiss[0])+
		   (AxePere->Meris->DirCroiss[1]*AxePere->Meris->DirCroiss[1]));
  if (NorVProjHor<Epsilon)
    {
      VAxeRot[0]=1.0; /* Vecteur initial vertical */
      VAxeRot[1]=0.0;
      VAxeRot[2]=0.0; /* Vecteur (1,0,0) choisi pour axe de rotation */
    }
  else
    {
      VAxeRot[0]=AxePere->Meris->DirCroiss[1]/NorVProjHor;
      VAxeRot[1]=-AxePere->Meris->DirCroiss[0]/NorVProjHor;
      VAxeRot[2]=0.0;
    }
  /* On fait tourner DirCroiss autour de VAxeRot d'un angle d'insertion */
  AngRot=TireAngI(TypeFils);
  RotVect(AngRot,VAxeRot,AxePere->Meris->DirCroiss,RotDirCroiss);

  /* On fait tourner RotDirCroiss autour de DirCroiss d'un angle radial */
  AngRot=TireAngRad();
  RotVect(AngRot,AxePere->Meris->DirCroiss,RotDirCroiss,DirFils);
} /* Fonction OrienteRamif */
/****************************************************************************/
/****************************************************************************/
void RamifieAxe(PTAxe AxePere)
{
  PTAxe NouvAxe;
  int TypeRamif;
  r3 OrigRamif, DirRamif;

  /* Décremente la distance au dernier primordium initié */
  AxePere->Meris->DistPrimInit=(AxePere->Meris->DistPrimInit)
    -DistInterRamifMeris(AxePere->Meris,Sol);

  /* Calcul de attributs d'une ramification */
  TypeRamif=TireTypeFils(AxePere);             /* Le type de son méristème */
  OrigineRamif(AxePere,OrigRamif);             /* Sa position */
  OrienteRamif(AxePere,TypeRamif,DirRamif);    /* Sa direction */

  NouvAxe=InitialiseAxe(SR->NbAxeForm+1,TypeRamif,OrigRamif,DirRamif,AxePere);

  AjouteAxeSR(SR,NouvAxe);

  AjouteNoeudLatAxe(AxePere,NouvAxe->PremNoeud);

} /* Fonction RamifieAxe */
/****************************************************************************/
/****************************************************************************/
void OrigineReit(PTAxe AxePere, r3 OrigineFils)
{
  /* Calcule le point d'origine d'une reiteration */
  OrigineFils[0]=AxePere->Meris->Coord[0];
  OrigineFils[1]=AxePere->Meris->Coord[1];
  OrigineFils[2]=AxePere->Meris->Coord[2];
} /* Fonction OrigineReit */
/****************************************************************************/
/****************************************************************************/
void OrienteReit(PTAxe AxePere, int TypeFils, r3 DirFils)
{
  /* Calcule la direction d'une reiteration */
  r3 VAxeRot,RotDirCroiss;
  double NorVProjHor,AngRot;

  /* Calcul de la norme de la projection direction sur plan horizontal */
  NorVProjHor=sqrt((AxePere->Meris->DirCroiss[0]*AxePere->Meris->DirCroiss[0])+
		   (AxePere->Meris->DirCroiss[1]*AxePere->Meris->DirCroiss[1]));
  if (NorVProjHor<Epsilon)
    {
      VAxeRot[0]=1.0; /* Vecteur initial vertical */
      VAxeRot[1]=0.0;
      VAxeRot[2]=0.0; /* Vecteur (1,0,0) choisi pour axe de rotation */
    }
  else
    {
      VAxeRot[0]=AxePere->Meris->DirCroiss[1]/NorVProjHor;
      VAxeRot[1]=-AxePere->Meris->DirCroiss[0]/NorVProjHor;
      VAxeRot[2]=0.0;
    }
  /* On fait tourner DirCroiss autour de VAxeRot d'un angle d'insertion */
  AngRot=TireAngIReit(TypeFils);
  RotVect(AngRot,VAxeRot,AxePere->Meris->DirCroiss,RotDirCroiss);

  /* On fait tourner RotDirCroiss autour de DirCroiss d'un angle generatrice */
  AngRot=TireAngRad();
  RotVect(AngRot,AxePere->Meris->DirCroiss,RotDirCroiss,DirFils);
} /* Fonction OrienteReit */
/****************************************************************************/
/****************************************************************************/
void ReitereAxe(PTAxe AxePere)
{
  PTAxe NouvAxe;
  int NumReit, NbReit;
  r3 Origine, DirReit;

  if ((AxePere->Meris->Mature)&&(!AxePere->Meris->Senile)&&
      (FRandUnif()<P_ProbReit[AxePere->Meris->Type]))
    {
      NbReit=IRandUnif(P_NbReitMax[AxePere->Meris->Type]-P_NbReitMin[AxePere->Meris->Type])
	+P_NbReitMin[AxePere->Meris->Type];
      for (NumReit=1; NumReit<=NbReit; NumReit++)
	{
	  OrigineReit(AxePere,Origine);
	  OrienteReit(AxePere,AxePere->Meris->Type,DirReit);
	  NouvAxe=InitialiseAxe(SR->NbAxeForm+1,AxePere->Meris->Type,Origine,DirReit,AxePere);
	  AjouteAxeSR(SR,NouvAxe);
	  AjouteNoeudLatAxe(AxePere,NouvAxe->PremNoeud);
	}
      AxePere->Meris->DistPrimInit=0.0;
      SenesceMeris(AxePere->Meris);
    }

} /* Fonction ReitereAxe */
/****************************************************************************/
/****************************************************************************/
void OrigineEmission(PTAxe NouvAxe)
{
  NouvAxe->Meris->Coord[0]=SR->Origine[0];
  NouvAxe->Meris->Coord[1]=SR->Origine[1];
  NouvAxe->Meris->Coord[2]=SR->Origine[2];
} /* Fonction OrigineEmission */
/****************************************************************************/
/****************************************************************************/
void OrienteEmission(PTAxe NouvAxe, int Num)
{
  double AngRot,AngI;
  r3 VInit;

  AngI=TireAngI(NouvAxe->Meris->Type);
  VInit[0]=sin(AngI);
  VInit[1]=0.0;
  VInit[2]=cos(AngI);

  AngRot=SR->AngDep+(2*Pi*Num/P_NbAxeOG);
  RotZ(VInit,NouvAxe->Meris->DirCroiss,AngRot);
} /* Fonction OrienteEmission */
/****************************************************************************/
/****************************************************************************/
void DeveloppeSR(PTSysRac SR)
{
  /* Développement, croissance, ramification et réitération de chaque axe */
  PTAxe AxeCour;

  AxeCour=SR->PremAxe;
  while (AxeCour!=NULL)
    {
      DeveloppeAxe(AxeCour);
      AllongeAxe(AxeCour);
      while (AxeRamifiable(AxeCour)) RamifieAxe(AxeCour);
      if (SR->ReitPossible) ReitereAxe(AxeCour);
      AxeCour=AxeCour->Suivant;
    }

/*  printf(" NbRac : %6d \n",SR->NbAxeForm);*/

}  /* Fonction DeveloppeSR */
/****************************************************************************/
/****************************************************************************/
void MortaliteSR(PTSysRac SR)
{
  PTAxe AxeCour, AxeAEnlever;
  int Necrose;

  /* Calcul de la mortalite sur l'ensemble des Axes */
  AxeCour=SR->DernAxe;
  while (AxeCour!=NULL)
    {
      if ((AxeCour->Meris->Senile)&&
	  (Temps-(DernNoeudAxe(AxeCour)->JourForm))>P_DureeNecrose[AxeCour->Meris->Type])
	{ /* L'axe est necrose */
	  Necrose=1;
	  AffecValNecroseAxe(AxeCour, Necrose);
	}
      else
	{  /* L'axe concerne n'est pas necrose */
	  Necrose=0;
	  AffecValNecroseAxe(AxeCour, Necrose);
	  /* Et les noeuds en amont ne sont pas necroses non plus */
	  AffecValNecroseAmont(AxeCour, Necrose);
	}

      AxeCour=AxeCour->Precedent;
    }

  /* Deuxieme passage pour re-specifier l'ensemble non necrose */
  AxeCour=SR->DernAxe;
  while (AxeCour!=NULL)
    {
      if ((AxeCour->Meris->Senile)&&
	  (Temps-(DernNoeudAxe(AxeCour)->JourForm))>P_DureeNecrose[AxeCour->Meris->Type])
	{ /* L'axe est necrose */
	}
      else
	{  /* L'axe concerne n'est pas necrose */
	  Necrose=0;
	  AffecValNecroseAxe(AxeCour, Necrose);
	  /* Et les noeuds en amont ne sont pas necroses non plus */
	  AffecValNecroseAmont(AxeCour, Necrose);
	}

      AxeCour=AxeCour->Precedent;
    }

  /* Calcul de l'elagage, enlevement des axes necroses */
  AxeCour=SR->DernAxe;
  while (AxeCour!=NULL)
    {
      if (AxeToutNecrose(AxeCour))
	{
	  AxeAEnlever=AxeCour;
	  AxeCour=AxeCour->Precedent;
	  EnleveAxeSR(SR,AxeAEnlever);
	}
      else AxeCour=AxeCour->Precedent;
    }

}  /* Fonction MortaliteSR */
/****************************************************************************/
/****************************************************************************/
void CroissanceRadialeSR(PTSysRac SR, float CoeffCroiss)
{
  PTAxe AxeCour;
  float Diam;

  /* Premier passage, initialisation aux diametres primaires */
  AxeCour=SR->DernAxe;
  while (AxeCour!=NULL)
    {
      Diam=P_DiamPrim[AxeCour->Meris->Type];
      AffecValDiamAxe(AxeCour, Diam);
      AxeCour=AxeCour->Precedent;
    }

  /* Deuxieme passage, avec increment des diametres */
  AxeCour=SR->DernAxe;
  while (AxeCour!=NULL)
    {
      /* les noeuds en amont sont incrementes */
      Diam=P_DiamPrim[AxeCour->Meris->Type];
      IncreValDiamAmont(AxeCour, Diam, CoeffCroiss);
      AxeCour=AxeCour->Precedent;
    }

}  /* Fonction CroissanceRadialeSR */
/****************************************************************************/
/****************************************************************************/
void ImprimeNd(PTNoeud Nd, long int NumAxe)
{
  long int SuivSF,SuivSP,Pere;

  if (Nd->SuivSFils==NULL) SuivSF=0;
  else SuivSF=Nd->SuivSFils->Num;

  if (Nd->SuivSPere==NULL) SuivSP=0;
  else SuivSP=Nd->SuivSPere->Num;

  if (Nd->Pere==NULL) Pere=-9;
  else Pere=Nd->Pere->Num;


  fprintf(FNoeud,"%5li %5i %2i %2i %5li %5li %5li %5li %7.2f %7.2f %7.2f %7.2f\n",
	  Nd->Num,Nd->JourForm,Nd->Fils->Meris->TypeDest,Nd->Necrose,NumAxe,SuivSF,SuivSP,Pere,
	  Nd->Diametre,Nd->Pos[0],Nd->Pos[1],Nd->Pos[2]);
}  /* Fonction ImprimeNd */
/****************************************************************************/
/****************************************************************************/
void ImprimeEnteteSR(void)
{   /* Imprime l'entête du fichier contenant le système racinaire */

  fprintf(FNoeud,"NumNd Jour Type Nec NumAxe SuivSF SuivSP Pere Diam     X       Y       Z\n");

}  /* Fonction ImprimeEnteteSR */
/****************************************************************************/
/****************************************************************************/
void ImprimeSansRamifAxe(PTAxe Axe)
{   /* Imprime les noeuds d'un axe sans y inclure les noeuds de branchement */

  PTNoeud NdCour;

  NdCour=Axe->PremNoeud;
  if (NdCour->SuivSFils!=NULL)
    {
      ImprimeNd(NdCour,Axe->Num); /* Ecriture du premier noeud dès lors qu'il a un suivant */
      NdCour=NdCour->SuivSFils;
      while ((NdCour->SuivSFils!=NULL)||(NdCour->SuivSPere!=NULL))
	{
	  if (NdCour->SuivSPere!=NULL) { NdCour=NdCour->SuivSPere; } /* On passe, c'est ramif */
	  else
	    {
	      ImprimeNd(NdCour,Axe->Num); /* Ecriture des noeuds intermédiaires non branchés */
	      NdCour=NdCour->SuivSFils;
	    }
	}
      ImprimeNd(NdCour,Axe->Num); /* Ecriture du dernier noeud */
    }

}  /* Fonction ImprimeSansRamifAxe */
/****************************************************************************/
/****************************************************************************/
void ImprimeAvecRamifAxe(PTAxe Axe)
{   /* Imprime les noeuds d'un axe en incluant les noeuds de branchement */

  PTNoeud NdCour;

  NdCour=Axe->PremNoeud;
  ImprimeNd(NdCour,Axe->Num);
  if (NdCour->SuivSFils!=NULL)
    {
      NdCour=NdCour->SuivSFils;
      ImprimeNd(NdCour,Axe->Num);
      while ((NdCour->SuivSFils!=NULL)||(NdCour->SuivSPere!=NULL))
	{
	  if (NdCour->SuivSPere!=NULL) NdCour=NdCour->SuivSPere;
	  else NdCour=NdCour->SuivSFils;
	  ImprimeNd(NdCour,Axe->Num);
	}
    }
}  /* Fonction ImprimeAvecRamifAxe */
/****************************************************************************/
/****************************************************************************/
void ImprimeAvecRamifSR(PTSysRac SR)
{   /* Imprime le système racinaire avec les noeuds de ramification */

  PTAxe AxeCour;

  ImprimeEnteteSR();

  AxeCour=SR->PremAxe;
  while (AxeCour!=NULL)
    {
      ImprimeAvecRamifAxe(AxeCour);
      AxeCour=AxeCour->Suivant;
    }
}  /* Fonction ImprimeAvecRamifSR */
/****************************************************************************/
/****************************************************************************/
void ImprimeSansRamifSR(PTSysRac SR)
{   /* Imprime le système racinaire sansles noeuds de ramification */

  PTAxe AxeCour;

  ImprimeEnteteSR();

  AxeCour=SR->PremAxe;
  while (AxeCour!=NULL)
    {
      ImprimeSansRamifAxe(AxeCour);
      AxeCour=AxeCour->Suivant;
    }
}  /* Fonction ImprimeSansRamifSR */
/****************************************************************************/
/****************************************************************************/
void FermeFichiers(void)
{
  fclose(FNoeud);
  fclose(FPar);
  fclose(FSol);
}  /* Fonction FermeFichiers */
/****************************************************************************/

void init_roottype1_(r3 Orig)
{
  srand(time(NULL)); /* Initialisation du générateur aléatoire Mathieu*/
  OuvreFichiers();
  LitParam();
  LitSol();
  LitLim();
  
  SR=InitialiseSR(Orig);
  InstalleSR(SR);
}


void iterate_roottype1_(int* p_duree)
{
  P_Duree=*p_duree;
  while (Temps<P_Duree) {
    Temps=Temps+DeltaT;
   printf("Temps : %3d %3d\n",Temps, P_Duree);
    
    /* Calcule l'état du système racinaire en terme de réitération */
    EtatReiterationSR(SR);
    
    /* Développement du système racinaire */
    printf("DeveloppeSR\n");
    DeveloppeSR(SR);
    
    /* Croissance radiale du système racinaire */
printf("CroissanceRadialeSR\n");
    CroissanceRadialeSR(SR, P_CoeffCroissRad);
    
    /* Mortalite du système racinaire */
printf("MortaliteSR\n");
    MortaliteSR(SR);
  }
}


void finish_roottype1_()
{
ImprimeAvecRamifSR(SR);
  FermeFichiers();
}

