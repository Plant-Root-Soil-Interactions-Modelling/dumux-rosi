#include <stdlib.h>
#include <stdio.h>
#include "RootTyp1boite.h"
#include "RootTyp1_interface.h"

static int *new_num=0;
static PTNoeud *node_ptr=0;
static int node_no=0;


void number_of_nodes_(int *n_nodes, int *n_axes_crea, int *n_axes_del)
{
  *n_axes_crea=SR->NbAxeForm;
  *n_axes_del=SR->NbAxeSup;
  *n_nodes=node_no=SR->NbNoeudForm;
}


void extract_nodes_(int *node, int *prev_node, int *ordn, double *xn, 
                    double *yn, double *zn, double *diamn, int *agen, int *axeF,
		      int *orda, double *xa, double *ya, double *za, double *diama, 
		      double *agea, int *prev_nodea, int *axenum, int *axenf,
		      int *node_corr)
{
  PTAxe p_axe;
  PTNoeud p_nodea;
  PTNoeud p_node;
  PTMeristeme p_meris;
  int no_first_nodes=0;
  int no_axe=0;
  int totaxe=0;
  int i, n, num;


  /* extract unique nodes */
  if(node_ptr) {
    free(node_ptr);
  }
  node_ptr=(PTNoeud*) malloc( (SR->NbNoeudForm+1) * sizeof(PTNoeud));
  node_ptr[0]=0;
  new_num=(int*) malloc( (SR->NbNoeudForm+1) * sizeof(int));
  new_num[0]=0;
  *node_corr=node_no;
  p_axe=SR->PremAxe;
  n=0;
  while(p_axe) {  /* loop over all axes */
      p_nodea=p_axe->DernNoeud;
      p_node=p_axe->PremNoeud;
      p_meris=p_axe->Meris;

    /* apex information*/
	
      if(p_nodea->Num!=p_node->Num && p_nodea->Prec->Num!=p_node->Num) {
          totaxe++;
          axeF[totaxe-1]=totaxe;
          xa[totaxe-1]=p_meris->Coord[0];
          ya[totaxe-1]=p_meris->Coord[1];
          za[totaxe-1]=p_meris->Coord[2];
          orda[totaxe-1]=p_meris->TypeDest;/* new rt version */
          diama[totaxe-1]=p_meris->Diametre;
          agea[totaxe-1]=p_meris->Age;

          if (totaxe>1)
              p_node=p_node->SuivSFils;
          while(p_node) {  /* loop over one axe */
	       if(p_node->Prec) {
                  prev_node[n]=new_num[p_node->Prec->Num];}	/*Couvreur may 2010*/
	       else /* close to the seed*/
	           prev_node[n]=0;
              new_num[p_node->Num]=n+1;
	       node[n]=new_num[p_node->Num];
	       axenf[n]=totaxe;
              xn[n]=p_node->Pos[0];
              yn[n]=p_node->Pos[1];
              zn[n]=p_node->Pos[2];
	       diamn[n]=p_node->Diametre;
	       ordn[n]=p_meris->TypeDest;/* new rt version */
              agen[n]=p_node->JourForm;
	       n++;

	       if(n>node_no) error(1);
              node_ptr[n]=p_node;
              if(p_nodea->Num==p_node->Num) {       
	           break; /* leave axe loop */
              }
              p_node=p_node->SuivSAxe;
          }  /* end of loop over one axe */
          prev_nodea[totaxe-1]=new_num[p_nodea->Num];		/*Couvreur may 2010*/
      }
      p_axe=p_axe->Suivant;
  }  /* end of loop over all axes */
  *axenum=totaxe;
  if(n!=node_no) {
      printf("%d nodes have been deleted",node_no-n);
      *node_corr=n;
      node_no=n;
  }
}


void write_node_data_(double *wc)
{
  int n;
  for(n=0;n<node_no;n++) {
    node_ptr[n+1]->wc=wc[n];
  }
}


void read_node_data_(double *wc)
{
  int n;
  for(n=0;n<node_no;n++) {
    wc[n]=node_ptr[n+1]->wc;
  }
}


void close_c_interface_()
{
  if(node_ptr) {
    free(node_ptr);    
    node_ptr=0;
  }
}


void error(int type)
{
  switch(type) {
  case 1:
    printf("Error: Invalid number of nodes in extract_nodes().\n");
    break;
  case 2:
    printf("Error: Number of nodes in extract_nodes() does not match.\n");
    break;
  }
  close_c_interface_();
  exit(type);
}
