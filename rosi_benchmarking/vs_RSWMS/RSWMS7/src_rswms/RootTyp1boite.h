#ifndef ROOTTYP1
#define ROOTTYP1

#define NBTYPMAX 13  /* Nbre maximal de types */
#define NBHORMAX 30 /* Nombre d'horizons de sol */
#define NBVAGMAX 12  /* Nombre maximal de vagues de réitération */
#define MAXLINE 120  /* Longueur maxi de la ligne dans fichiers texte */

typedef double r2[2];  /* Tableau 2D */
typedef double r3[3];  /* Tableau 3D */

typedef struct sysrac *PTSysRac;
typedef struct axe *PTAxe;
typedef struct meristeme *PTMeristeme;
typedef struct noeud *PTNoeud;

typedef struct sysrac /* Ensemble d'axes */
{
  long int NbAxeForm;  /* Nombre d'axes formés */
  long int NbAxeSup;   /* Nombre d'axes supprimés */
  long int NbNoeudForm;    /* Nombre de noeuds formés */
  int NumVReitCourante; /* Numéro de vague de réitération courante */
  int ReitPossible;     /* Réitération possible : 1 ou 0 */
  double AngDep; /* Orientation */
  r3 Origine; /* Position de l'origine */
  PTAxe PremAxe;       /* Premier axe du système (accès à la liste) */
  PTAxe DernAxe;       /* Dernier axe produit */
} TSysRac ;

typedef struct meristeme /* Méristème apical, ou pointe de chaque racine */
{
  double DistPrimInit; /* Distance de l'apex au dernier primordium initié */
  r3 Coord;           /* Coordonnées de l'apex */
  r3 DirCroiss;       /* Direction de croissance */
  r3 DirInit;         /* Direction initiale */
  double Age;          /* Age du méristème */
  double Diametre;     /* Diamètre au niveau de la pointe */
  r2 PCroiss;         /* Caractéristiques de croissance potentielle du méristème */
  double PRamif;       /* Distance inter-ramif potentielle */
  int Type;           /* Type (état morphogénétique) de méristème : de 0 a 7 */
  int TypeDest;       /* Type à faire apparaître dans le fichier noeud de sortie */
  int Senile;         /* Sénile ?, ou actif ... */
  int Mature;         /* Mature ?, ou primordium ... */
} TMeristeme ;

typedef struct axe /* Ensemble constitué d'un méristème et liste de noeuds */
{
  long int Num;      /* Numéro de l'axe */
  int NbNoeud;       /* Nombre de noeuds */
  PTMeristeme Meris; /* Méristème apical */
  PTAxe Suivant;     /* Suivant de la liste */
  PTAxe Precedent;   /* Precedent de la liste */
  PTAxe Pere;        /* Axe père, sur lequel celui-ci est branché */
  PTNoeud PremNoeud; /* Premier noeud de l'axe, sa base */
  PTNoeud DernNoeud; /* Dernier noeud de l'axe, apical */
} TAxe ;

typedef struct noeud
{
  long int Num;      /* Numéro d'ordre de création */
  int JourForm;      /* Date de formation (en jours) */
  double Diametre;   /* Diametre au niveau du noeud */
  double wc;         /* water content of the node (HH) */ 
  r3 Pos;            /* Position dans l'espace */
  PTNoeud SuivSPere; /* Suivant sur le père (quand branchement) */
  PTNoeud SuivSFils; /* Suivant sur le fils, axe d'appartenance */
  PTNoeud Prec;      /* Precedent */
  PTNoeud SuivSAxe;  /* Suivant sur l'axe */
  PTAxe Pere;        /* Axe père, pointe sur NULL quand non branchement */
  PTAxe Fils;        /* Axe fils, d'appartenance du noeud */
  int Necrose;       /* Necrose ? 0 : non; 1 : oui */
} TNoeud ;

typedef struct horizon
{
  float Croiss;  /* Coefficient de croissance, compris entre 0 et 1 */
  float Ramif;   /* Coefficient multiplicateur de distance inter-ramif  >1 */
  float ICMeca;  /* Intensité de la contrainte mécanique */
  int OCMeca;    /* Orientation de la contrainte mécanique (O iso, ou 1 vert) */
} THorizon ;

typedef THorizon TSol[NBHORMAX];  /* Sol, tableau d'horizons */

/* Variables globales diverses */
extern int Temps;  /* Le temps, en jours */
extern r3 Orig;    /* Position d'origine de l'organe générateur */

extern PTSysRac SR;  /* Le système racinaire */
extern TSol Sol;     /* Le sol */
#endif
