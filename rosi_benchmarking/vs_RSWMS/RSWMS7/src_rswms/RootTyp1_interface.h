#ifndef ROOTTYP1_INTERFACE
#define ROOTTYP1_INTERFACE

void error(int type);
void number_of_nodes_(int *n_nodes, int *n_axes_crea, int *n_axes_del);
void extract_nodes_(int *node, int *prev_node, int *ordn, double *xn, double *yn, double *zn,
 double *diamn, int *agen, int *axeF, int *orda, double *xa, double *ya, double *za, double *diama, 
 double *agea, int *prev_nodea, int *totaxe, int *axenf, int *node_no);
void write_node_data_(double *wc);
void read_node_data_(double *wc);
void close_c_interface_();
#endif
