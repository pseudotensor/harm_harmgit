
/*! \file coord.funcdeclare.h
    \brief Coordinates stuff (both user and general)    
*/


extern void set_points(void);
// coordinate stuff
extern void set_coord_parms(int defcoordlocal);
extern void set_coord_parms_nodeps(int defcoordlocal);
extern void set_coord_parms_deps(int defcoordlocal);
extern void write_coord_parms(int defcoordlocal);
extern void read_coord_parms(int defcoordlocal);
extern void coord(int i, int j, int k, int loc, FTYPE *X);
extern void coord_ijk(int i, int j, int k, int loc, FTYPE *X);
extern void coord_free(int i, int j, int k, int loc, FTYPE *X);
extern void coordf(FTYPE i, FTYPE j, FTYPE k, int loc, FTYPE *X);

extern void bl_coord(FTYPE *X, FTYPE *V);
extern void bl_coord_ijk(int i, int j, int k, int loc, FTYPE *V);
extern void bl_coord_ijk_2(int i, int j, int k, int loc, FTYPE *X, FTYPE *V);
extern void bl_coord_coord(int i, int j, int k, int loc, FTYPE *X, FTYPE *V);

extern void dxdxprim(FTYPE *X, FTYPE *V, FTYPE (*dxdxp)[NDIM]);
extern void dxdxprim_ijk(int i, int j, int k, int loc, FTYPE (*dxdxp)[NDIM]);
extern void dxdxprim_ijk_2(struct of_geom *ptrgeom, FTYPE *X, FTYPE *V, FTYPE (*dxdxp)[NDIM]);

extern void idxdxprim(FTYPE (*dxdxp)[NDIM], FTYPE (*idxdxp)[NDIM]);
extern void idxdxprim_ijk(int i, int j, int k, int loc, FTYPE (*idxdxp)[NDIM]);
extern void idxdxprim_ijk_2(struct of_geom *ptrgeom, FTYPE *X, FTYPE *V, FTYPE (*idxdxp)[NDIM]);


extern int setihor(void);
extern FTYPE setRin(int ihor);

extern int is_inside_surface(int dir, int ii, int jj, int kk, int pp);
extern int is_on_surface(int dir, int ii, int jj, int kk, int pp);
