
/*! \file global.fieldmacros.h
    \brief Magnetic field and divB=0 related macros
    
*/



#if(NEWMETRICSTORAGE==1)

#if(MCOORD!=CARTMINKMETRIC)
/// e.g. MYGDET(i,j,k,CENT)
/// GODMARK: assumes field EOMs have WHICHEOM==WITHGDET or NOGDETB?=0
#define MYGDET(i,j,k,p) (GLOBALMETMACP1A0(compgeom,p,i,j,k).gdet)
#else
#define MYGDET(i,j,k,p) (GLOBALMETMACP1A0(compgeom,p,0,0,0).gdet)
#endif


#else // old method

#if(MCOORD!=CARTMINKMETRIC)
/// e.g. MYGDET(i,j,k,CENT)
#define MYGDET(i,j,k,p) (GLOBALMETMACP1A0(gdet,p,i,j,k_))
#else
#define MYGDET(i,j,k,p) (GLOBALMETMACP1A0(gdet,p,0,0,0))
#endif


#endif // end if old method


//////////////////////////////////
///
/// below 4 macros used for vpot2field() in initbase.c
#define FgCORN(F,i,j,k) (MAC(F,i,j,k))
/// AVG_x functions that average in x direction
/// typically operates on CORN1,2,3 quantities to get face quantities
#define AVGCORN_1(F,i,j,k) (0.5*(FgCORN(F,ip1mac(i),j,k) + FgCORN(F,i,j,k) ))
#define AVGCORN_2(F,i,j,k) (0.5*(FgCORN(F,i,jp1mac(j),k) + FgCORN(F,i,j,k) ))
#define AVGCORN_3(F,i,j,k) (0.5*(FgCORN(F,i,j,kp1mac(k)) + FgCORN(F,i,j,k) ))

/// do-nothing functions
#define NOAVGCORN_1(F,i,j,k) (FgCORN(F,i,j,k))
#define NOAVGCORN_2(F,i,j,k) (FgCORN(F,i,j,k))
#define NOAVGCORN_3(F,i,j,k) (FgCORN(F,i,j,k))


////////////////////////////////
///
/// below macros are similar to those farther below for divb, but no gdet stuff
#define FgN(F,i,j,k,pl) (MACP0A1(F,i,j,k,pl))

/// AVGN_x functions that average in x direction
/// typically operates on centered quantities to get face quantities OR on face quantities to get edge quantities
#define AVGN_1(F,i,j,k,pl) (0.5*(FgN(F,i,j,k,pl) + FgN(F,im1mac(i),j,k,pl) ))
#define AVGN_2(F,i,j,k,pl) (0.5*(FgN(F,i,j,k,pl) + FgN(F,i,jm1mac(j),k,pl) ))
#define AVGN_3(F,i,j,k,pl) (0.5*(FgN(F,i,j,k,pl) + FgN(F,i,j,km1mac(k),pl) ))

#define NOAVGN_1(F,i,j,k,pl) (FgN(F,i,j,k,pl))
#define NOAVGN_2(F,i,j,k,pl) (FgN(F,i,j,k,pl))
#define NOAVGN_3(F,i,j,k,pl) (FgN(F,i,j,k,pl))

/// just average of 4 quantities that eventually lives on (0.5,0,0) (i.e. CORN1)
#define AVGN_for1(F,i,j,k,pl) (0.5*(AVGN_2(F,i,j,k,pl)+AVGN_2(F,i,j,km1mac(k),pl))) // 2 and 3 directions averaged
/// lives on (0,0.5,0) (i.e. CORN2)
#define AVGN_for2(F,i,j,k,pl) (0.5*(AVGN_1(F,i,j,k,pl)+AVGN_1(F,i,j,km1mac(k),pl))) // 1 and 3 directions averaged
/// lives on (0,0,0.5) (i.e. CORN3)
#define AVGN_for3(F,i,j,k,pl) (0.5*(AVGN_2(F,i,j,k,pl)+AVGN_2(F,im1mac(i),j,k,pl))) // 1 and 2 directions averaged

#define NOAVGN_for1(F,i,j,k,pl) (NOAVGN_1(F,i,j,k,pl))
#define NOAVGN_for2(F,i,j,k,pl) (NOAVGN_2(F,i,j,k,pl))
#define NOAVGN_for3(F,i,j,k,pl) (NOAVGN_3(F,i,j,k,pl))


////////////////////////////////
///
/// below macros are similar to those farther below for divb, but no gdet stuff
#define absFgN(F,i,j,k,pl) (fabs(MACP0A1(F,i,j,k,pl)))

/// AVGN_x functions that average in x direction
/// typically operates on centered quantities to get face quantities OR on face quantities to get edge quantities
#define absAVGN_1(F,i,j,k,pl) (0.5*(absFgN(F,i,j,k,pl) + absFgN(F,im1mac(i),j,k,pl) ))
#define absAVGN_2(F,i,j,k,pl) (0.5*(absFgN(F,i,j,k,pl) + absFgN(F,i,jm1mac(j),k,pl) ))
#define absAVGN_3(F,i,j,k,pl) (0.5*(absFgN(F,i,j,k,pl) + absFgN(F,i,j,km1mac(k),pl) ))

#define absNOAVGN_1(F,i,j,k,pl) (absFgN(F,i,j,k,pl))
#define absNOAVGN_2(F,i,j,k,pl) (absFgN(F,i,j,k,pl))
#define absNOAVGN_3(F,i,j,k,pl) (absFgN(F,i,j,k,pl))

/// just average of 4 quantities that eventually lives on (0.5,0,0) (i.e. CORN1)
#define absAVGN_for1(F,i,j,k,pl) (0.5*(absAVGN_2(F,i,j,k,pl)+absAVGN_2(F,i,j,km1mac(k),pl))) // 2 and 3 directions averaged
/// lives on (0,0.5,0) (i.e. CORN2)
#define absAVGN_for2(F,i,j,k,pl) (0.5*(absAVGN_1(F,i,j,k,pl)+absAVGN_1(F,i,j,km1mac(k),pl))) // 1 and 3 directions averaged
/// lives on (0,0,0.5) (i.e. CORN3)
#define absAVGN_for3(F,i,j,k,pl) (0.5*(absAVGN_2(F,i,j,k,pl)+absAVGN_2(F,im1mac(i),j,k,pl))) // 1 and 2 directions averaged

#define absNOAVGN_for1(F,i,j,k,pl) (absNOAVGN_1(F,i,j,k,pl))
#define absNOAVGN_for2(F,i,j,k,pl) (absNOAVGN_2(F,i,j,k,pl))
#define absNOAVGN_for3(F,i,j,k,pl) (absNOAVGN_3(F,i,j,k,pl))


///////////////////////////////////////////
///
/// below macros used for divb definition (with gdet)
#define Fg(F,i,j,k,pl) (MACP0A1(F,i,j,k,pl)*MYGDET(i,j,k,CENT))
#define Fgface(F,i,j,k,pl) (MACP0A1(F,i,j,k,pl)*MYGDET(i,j,k,FACE1+pl-B1)) // assumes B1,B2,B3 in order for pl GODMARK

/// AVG_x functions that average in x direction (includes geometry)
/// typically operates on centered quantities to get face quantities OR on face quantities to get edge quantities
#define AVG_1(F,i,j,k,pl) (0.5*(Fg(F,i,j,k,pl) + Fg(F,im1mac(i),j,k,pl) ))
#define AVG_2(F,i,j,k,pl) (0.5*(Fg(F,i,j,k,pl) + Fg(F,i,jm1mac(j),k,pl) ))
#define AVG_3(F,i,j,k,pl) (0.5*(Fg(F,i,j,k,pl) + Fg(F,i,j,km1mac(k),pl) ))

#define NOAVG_1(F,i,j,k,pl) (Fg(F,i,j,k,pl))
#define NOAVG_2(F,i,j,k,pl) (Fg(F,i,j,k,pl))
#define NOAVG_3(F,i,j,k,pl) (Fg(F,i,j,k,pl))


#define NOAVGFACE_1(F,i,j,k,pl) (Fgface(F,i,j,k,pl))
#define NOAVGFACE_2(F,i,j,k,pl) (Fgface(F,i,j,k,pl))
#define NOAVGFACE_3(F,i,j,k,pl) (Fgface(F,i,j,k,pl))

/// just average of 4 quantities that eventually lives on (0.5,0,0) (i.e. CORN1)
#define AVG_for1(F,i,j,k,pl) (0.5*(AVG_2(F,i,j,k,pl)+AVG_2(F,i,j,km1mac(k),pl))) // 2 and 3 directions averaged
/// lives on (0,0.5,0) (i.e. CORN2)
#define AVG_for2(F,i,j,k,pl) (0.5*(AVG_1(F,i,j,k,pl)+AVG_1(F,i,j,km1mac(k),pl))) // 1 and 3 directions averaged
/// lives on (0,0,0.5) (i.e. CORN3)
#define AVG_for3(F,i,j,k,pl) (0.5*(AVG_2(F,i,j,k,pl)+AVG_2(F,im1mac(i),j,k,pl))) // 1 and 2 directions averaged

#define NOAVG_for1(F,i,j,k,pl) (NOAVG_1(F,i,j,k,pl))
#define NOAVG_for2(F,i,j,k,pl) (NOAVG_2(F,i,j,k,pl))
#define NOAVG_for3(F,i,j,k,pl) (NOAVG_3(F,i,j,k,pl))

#define NOAVGFACE_for1(F,i,j,k,pl) (NOAVGFACE_1(F,i,j,k,pl))
#define NOAVGFACE_for2(F,i,j,k,pl) (NOAVGFACE_2(F,i,j,k,pl))
#define NOAVGFACE_for3(F,i,j,k,pl) (NOAVGFACE_3(F,i,j,k,pl))


////////////////////////////////////////////
///
/// below macros used for divb definition (with gdet)
#define absFg(F,i,j,k,pl) (fabs(MACP0A1(F,i,j,k,pl)*MYGDET(i,j,k,CENT)))
#define absFgface(F,i,j,k,pl) (fabs(MACP0A1(F,i,j,k,pl)*MYGDET(i,j,k,FACE1+pl-B1))) // assumes B1,B2,B3 in order for pl GODMARK

/// absAVG_x functions that average in x direction
/// typically operates on centered quantities to get face quantities OR on face quantities to get edge quantities
#define absAVG_1(F,i,j,k,pl) (0.5*(absFg(F,i,j,k,pl) + absFg(F,im1mac(i),j,k,pl) ))
#define absAVG_2(F,i,j,k,pl) (0.5*(absFg(F,i,j,k,pl) + absFg(F,i,jm1mac(j),k,pl) ))
#define absAVG_3(F,i,j,k,pl) (0.5*(absFg(F,i,j,k,pl) + absFg(F,i,j,km1mac(k),pl) ))

#define absAVGFACE_1(F,i,j,k,pl) (0.5*(absFgface(F,i,j,k,pl) + absFgface(F,im1mac(i),j,k,pl) ))
#define absAVGFACE_2(F,i,j,k,pl) (0.5*(absFgface(F,i,j,k,pl) + absFgface(F,i,jm1mac(j),k,pl) ))
#define absAVGFACE_3(F,i,j,k,pl) (0.5*(absFgface(F,i,j,k,pl) + absFgface(F,i,j,km1mac(k),pl) ))

/// absAVG_x functions that average in x direction
/// typically operates on centered quantities to get face quantities OR on face quantities to get edge quantities
#define absNOAVG_1(F,i,j,k,pl) (absFg(F,i,j,k,pl))
#define absNOAVG_2(F,i,j,k,pl) (absFg(F,i,j,k,pl))
#define absNOAVG_3(F,i,j,k,pl) (absFg(F,i,j,k,pl))

#define absNOAVGFACE_1(F,i,j,k,pl) (absFgface(F,i,j,k,pl))
#define absNOAVGFACE_2(F,i,j,k,pl) (absFgface(F,i,j,k,pl))
#define absNOAVGFACE_3(F,i,j,k,pl) (absFgface(F,i,j,k,pl))

/// just average of 4 quantities that eventually lives on (0.5,0,0) (i.e. CORN1)
#define absAVG_for1(F,i,j,k,pl) (0.5*(absAVG_2(F,i,j,k,pl)+absAVG_2(F,i,j,km1mac(k),pl))) // 2 and 3 directions averaged
/// lives on (0,0.5,0) (i.e. CORN2)
#define absAVG_for2(F,i,j,k,pl) (0.5*(absAVG_1(F,i,j,k,pl)+absAVG_1(F,i,j,km1mac(k),pl))) // 1 and 3 directions averaged
/// lives on (0,0,0.5) (i.e. CORN3)
#define absAVG_for3(F,i,j,k,pl) (0.5*(absAVG_2(F,i,j,k,pl)+absAVG_2(F,im1mac(i),j,k,pl))) // 1 and 2 directions averaged

#define absNOAVG_for1(F,i,j,k,pl) (absAVG_1(F,i,j,k,pl))
#define absNOAVG_for2(F,i,j,k,pl) (absAVG_2(F,i,j,k,pl))
#define absNOAVG_for3(F,i,j,k,pl) (absAVG_3(F,i,j,k,pl))

#define absNOAVGFACE_for1(F,i,j,k,pl) (absAVGFACE_1(F,i,j,k,pl))
#define absNOAVGFACE_for2(F,i,j,k,pl) (absAVGFACE_2(F,i,j,k,pl))
#define absNOAVGFACE_for3(F,i,j,k,pl) (absAVGFACE_3(F,i,j,k,pl))






// #define DIVBCONDITION(p,i,j)
// if((i>=-1)&&(j>=-1)&&(startpos[2]+j!=0)&&(startpos[2]+j!=N2TOT))
//#define DIVBCONDITION(p,i,j) if((startpos[1]+i>0)&&(startpos[2]+j>0)&&(startpos[1]+i<totalsize[1])&&(startpos[2]+j<totalsize[2]))
//#define DIVBCONDITION(p,i,j) if((startpos[1]+i!=0)&&(startpos[1]+i!=totalsize[1])&&(startpos[2]+j!=0)&&(startpos[2]+j!=totalsize[2]))

// only needed for polar axis condition
//GODMARK
//#define DIVBCONDITION(p,i,j,k) if((startpos[2]+j!=0)&&(startpos[2]+j!=totalsize[2]))

// GODMARK: No, need for all boundaries
//#define DIVBCONDITION(p,i,j,k) if(((startpos[1]+i!=0)&&(startpos[1]+i!=totalsize[1]) || (!N1NOT1))&&((startpos[2]+j!=0)&&(startpos[2]+j!=totalsize[2]) || (!N2NOT1))&&((startpos[3]+k!=0)&&(startpos[3]+k!=totalsize[3]) || (!N3NOT1)))

/// only consider within bounds or ignore condition if no such dimension
#define DIVBCDIR1 ((i >= -N1BND+1 && i <= N1 - 1 + N1BND - 1) || (!N1NOT1))
#define DIVBCDIR2 ((j >= -N2BND+1 && j <= N2 - 1 + N2BND - 1) || (!N2NOT1))
#define DIVBCDIR3 ((k >= -N3BND+1 && k <= N3 - 1 + N3BND - 1) || (!N3NOT1))

#define DIVBCONDITION(p,i,j,k) if(DIVBCDIR1&&DIVBCDIR2&&DIVBCDIR3)

//#define DIVBCONDITION(p,i,j,k) if(1)




///////////////////////////////////
//
// FLUXCTTOTH DIVB
//
////////////////////////////////////
//(startpos[2]+j!=0)&&(startpos[2]+j!=totalsize[2])&&(startpos[3]+k!=0)&&(startpos[3]+k!=totalsize[3])

//////////////////////
/// PRMITIVE quantity
/// the below lines results in quantity at TRUE CORNER
#define DIVBDIFFFLUXCTTOTHPRIMx(p,i,j,k) (AVG_for1(p,i,j,k,B1)-AVG_for1(p,im1mac(i),j,k,B1))
#define DIVBDIFFFLUXCTTOTHPRIMy(p,i,j,k) (AVG_for2(p,i,j,k,B2)-AVG_for2(p,i,jm1mac(j),k,B2))
#define DIVBDIFFFLUXCTTOTHPRIMz(p,i,j,k) (AVG_for3(p,i,j,k,B3)-AVG_for3(p,i,j,km1mac(k),B3))

//#define DIVBNORMFLUXCTTOTHPRIMx(p,i,j,k) (AVG_for1(p,i,j,k,B1)+AVG_for1(p,im1mac(i),j,k,B1))
//#define DIVBNORMFLUXCTTOTHPRIMy(p,i,j,k) (AVG_for2(p,i,j,k,B2)+AVG_for2(p,i,jm1mac(j),k,B2))
//#define DIVBNORMFLUXCTTOTHPRIMz(p,i,j,k) (AVG_for3(p,i,j,k,B3)+AVG_for3(p,i,j,km1mac(k),B3))

#define DIVBNORMFLUXCTTOTHPRIMx(p,i,j,k) (absAVG_for1(p,i,j,k,B1)+absAVG_for1(p,im1mac(i),j,k,B1))
#define DIVBNORMFLUXCTTOTHPRIMy(p,i,j,k) (absAVG_for2(p,i,j,k,B2)+absAVG_for2(p,i,jm1mac(j),k,B2))
#define DIVBNORMFLUXCTTOTHPRIMz(p,i,j,k) (absAVG_for3(p,i,j,k,B3)+absAVG_for3(p,i,j,km1mac(k),B3))

// the below as dividing by too many \detg's since MYGDET in division and inside B
//#define DIVBNORMFLUXCTTOTHPRIM(p,i,j,k) (MAX(MAX(dx[1],dx[2]),dx[3])/(THIRD*MYGDET(i,j,k,CORNT)*fabs(DIVBNORMFLUXCTTOTHPRIMx(p,i,j,k)+DIVBNORMFLUXCTTOTHPRIMy(p,i,j,k)+DIVBNORMFLUXCTTOTHPRIMz(p,i,j,k)) +SMALL))
// removed MYGDET at corner.  This 1) removes issue with dividing by 0 and 2) is correct since DIVBNORMFLUXCTTOTHPRIM already has gdet
//#define DIVBNORMFLUXCTTOTHPRIM(p,i,j,k) (MAX(MAX(dx[1],dx[2]),dx[3])/(THIRD*fabs(DIVBNORMFLUXCTTOTHPRIMx(p,i,j,k)+DIVBNORMFLUXCTTOTHPRIMy(p,i,j,k)+DIVBNORMFLUXCTTOTHPRIMz(p,i,j,k)) +SMALL))
#define DIVBNORMFLUXCTTOTHPRIM(p,i,j,k) (1.0/(fabs(DIVBNORMFLUXCTTOTHPRIMx(p,i,j,k)/dx[1]+DIVBNORMFLUXCTTOTHPRIMy(p,i,j,k)/dx[2]+DIVBNORMFLUXCTTOTHPRIMz(p,i,j,k)/dx[3]) +SMALL))

#define DIVBFLUXCTTOTHPRIM(p,i,j,k) ((                                  \
                                      DIVBDIFFFLUXCTTOTHPRIMx(p,i,j,k)/dx[1] + DIVBDIFFFLUXCTTOTHPRIMy(p,i,j,k)/dx[2] + DIVBDIFFFLUXCTTOTHPRIMz(p,i,j,k)/dx[3] \
                                      )*(DIVBNORMFLUXCTTOTHPRIM(p,i,j,k)))



//////////////////////
/// CONSERVED quantity
/// the below lines results in quantity at TRUE CORNER
#define DIVBDIFFFLUXCTTOTHx(U,i,j,k) (AVGN_for1(U,i,j,k,B1)-AVGN_for1(U,im1mac(i),j,k,B1))
#define DIVBDIFFFLUXCTTOTHy(U,i,j,k) (AVGN_for2(U,i,j,k,B2)-AVGN_for2(U,i,jm1mac(j),k,B2))
#define DIVBDIFFFLUXCTTOTHz(U,i,j,k) (AVGN_for3(U,i,j,k,B3)-AVGN_for3(U,i,j,km1mac(k),B3))

//#define DIVBNORMFLUXCTTOTHx(U,i,j,k) (AVGN_for1(U,i,j,k,B1)+AVGN_for1(U,im1mac(i),j,k,B1))
//#define DIVBNORMFLUXCTTOTHy(U,i,j,k) (AVGN_for2(U,i,j,k,B2)+AVGN_for2(U,i,jm1mac(j),k,B2))
//#define DIVBNORMFLUXCTTOTHz(U,i,j,k) (AVGN_for3(U,i,j,k,B3)+AVGN_for3(U,i,j,km1mac(k),B3))

#define DIVBNORMFLUXCTTOTHx(U,i,j,k) (absAVGN_for1(U,i,j,k,B1)+absAVGN_for1(U,im1mac(i),j,k,B1))
#define DIVBNORMFLUXCTTOTHy(U,i,j,k) (absAVGN_for2(U,i,j,k,B2)+absAVGN_for2(U,i,jm1mac(j),k,B2))
#define DIVBNORMFLUXCTTOTHz(U,i,j,k) (absAVGN_for3(U,i,j,k,B3)+absAVGN_for3(U,i,j,km1mac(k),B3))

/// the below as dividing by too many \detg's since MYGDET in division and inside B
//#define DIVBNORMFLUXCTTOTH(U,i,j,k) (MAX(MAX(dx[1],dx[2]),dx[3])/(THIRD*MYGDET(i,j,k,CORNT)*fabs(DIVBNORMFLUXCTTOTHx(U,i,j,k)+DIVBNORMFLUXCTTOTHy(U,i,j,k)+DIVBNORMFLUXCTTOTHz(U,i,j,k)) +SMALL))
// removed MYGDET at corner.  This 1) removes issue with dividing by 0 and 2) is correct since DIVBNORMFLUXCTTOTH already has gdet
//#define DIVBNORMFLUXCTTOTH(U,i,j,k) (MAX(MAX(dx[1],dx[2]),dx[3])/(THIRD*fabs(DIVBNORMFLUXCTTOTHx(U,i,j,k)+DIVBNORMFLUXCTTOTHy(U,i,j,k)+DIVBNORMFLUXCTTOTHz(U,i,j,k)) +SMALL))
#define DIVBNORMFLUXCTTOTH(U,i,j,k) (1.0/(fabs(DIVBNORMFLUXCTTOTHx(U,i,j,k)/dx[1]+DIVBNORMFLUXCTTOTHy(U,i,j,k)/dx[2]+DIVBNORMFLUXCTTOTHz(U,i,j,k)/dx[3]) +SMALL))

#define DIVBFLUXCTTOTH(U,i,j,k) ((                                      \
                                  DIVBDIFFFLUXCTTOTHx(U,i,j,k)/dx[1] + DIVBDIFFFLUXCTTOTHy(U,i,j,k)/dx[2] + DIVBDIFFFLUXCTTOTHz(U,i,j,k)/dx[3] \
                                  )*(DIVBNORMFLUXCTTOTH(U,i,j,k)))




///////////////////////////////////
//
// FLUXCTSTAG DIVB (divb sits at center and assume p here is located at FACE)
//
///////////////////////////////////

///////////////////////
/// CONSERVED quantity
/// as input for divb calculation (otherwise have to change MYGDET so doesn't always use CENT
/// the below lines results in quantity at CENT
#define DIVBDIFFFLUXCTSTAGx(U,i,j,k) (NOAVGN_for1(U,ip1mac(i),j,k,B1)-NOAVGN_for1(U,i,j,k,B1))
#define DIVBDIFFFLUXCTSTAGy(U,i,j,k) (NOAVGN_for2(U,i,jp1mac(j),k,B2)-NOAVGN_for2(U,i,j,k,B2))
#define DIVBDIFFFLUXCTSTAGz(U,i,j,k) (NOAVGN_for3(U,i,j,kp1mac(k),B3)-NOAVGN_for3(U,i,j,k,B3))

#define DIVBNORMFLUXCTSTAGx(U,i,j,k) (absNOAVGN_for1(U,ip1mac(i),j,k,B1)+absNOAVGN_for1(U,i,j,k,B1))
#define DIVBNORMFLUXCTSTAGy(U,i,j,k) (absNOAVGN_for2(U,i,jp1mac(j),k,B2)+absNOAVGN_for2(U,i,j,k,B2))
#define DIVBNORMFLUXCTSTAGz(U,i,j,k) (absNOAVGN_for3(U,i,j,kp1mac(k),B3)+absNOAVGN_for3(U,i,j,k,B3))

//#define DIVBNORMFLUXCTSTAG(U,i,j,k) (MAX(MAX(dx[1],dx[2]),dx[3])/(THIRD*fabs(DIVBNORMFLUXCTSTAGx(U,i,j,k)+DIVBNORMFLUXCTSTAGy(U,i,j,k)+DIVBNORMFLUXCTSTAGz(U,i,j,k)) +SMALL))
#define DIVBNORMFLUXCTSTAG(U,i,j,k) (1.0/(fabs(DIVBNORMFLUXCTSTAGx(U,i,j,k)/dx[1]+DIVBNORMFLUXCTSTAGy(U,i,j,k)/dx[2]+DIVBNORMFLUXCTSTAGz(U,i,j,k)/dx[3]) +SMALL))

#define DIVBFLUXCTSTAG(U,i,j,k) ((                                      \
                                  DIVBDIFFFLUXCTSTAGx(U,i,j,k)/dx[1] + DIVBDIFFFLUXCTSTAGy(U,i,j,k)/dx[2] + DIVBDIFFFLUXCTSTAGz(U,i,j,k)/dx[3] \
                                  )*(DIVBNORMFLUXCTSTAG(U,i,j,k)))

//////////////////////
/// PRMITIVE quantity
/// as input for divb calculation (otherwise have to change MYGDET so doesn't always use CENT
/// the below lines results in quantity at CENT
#define DIVBDIFFFLUXCTSTAGPRIMx(p,i,j,k) (NOAVGFACE_for1(p,ip1mac(i),j,k,B1)-NOAVGFACE_for1(p,i,j,k,B1))
#define DIVBDIFFFLUXCTSTAGPRIMy(p,i,j,k) (NOAVGFACE_for2(p,i,jp1mac(j),k,B2)-NOAVGFACE_for2(p,i,j,k,B2))
#define DIVBDIFFFLUXCTSTAGPRIMz(p,i,j,k) (NOAVGFACE_for3(p,i,j,kp1mac(k),B3)-NOAVGFACE_for3(p,i,j,k,B3))

#define DIVBNORMFLUXCTSTAGPRIMx(p,i,j,k) (absNOAVGFACE_for1(p,ip1mac(i),j,k,B1)+absNOAVGFACE_for1(p,i,j,k,B1))
#define DIVBNORMFLUXCTSTAGPRIMy(p,i,j,k) (absNOAVGFACE_for2(p,i,jp1mac(j),k,B2)+absNOAVGFACE_for2(p,i,j,k,B2))
#define DIVBNORMFLUXCTSTAGPRIMz(p,i,j,k) (absNOAVGFACE_for3(p,i,j,kp1mac(k),B3)+absNOAVGFACE_for3(p,i,j,k,B3))

//#define DIVBNORMFLUXCTSTAGPRIM(p,i,j,k) (MAX(MAX(dx[1],dx[2]),dx[3])/(THIRD*fabs(DIVBNORMFLUXCTSTAGPRIMx(p,i,j,k)+DIVBNORMFLUXCTSTAGPRIMy(p,i,j,k)+DIVBNORMFLUXCTSTAGPRIMz(p,i,j,k)) +SMALL))
/// keep like-dimensional things together
#define DIVBNORMFLUXCTSTAGPRIM(p,i,j,k) (1.0/(fabs(DIVBNORMFLUXCTSTAGPRIMx(p,i,j,k)/dx[1]+DIVBNORMFLUXCTSTAGPRIMy(p,i,j,k)/dx[2]+DIVBNORMFLUXCTSTAGPRIMz(p,i,j,k)/dx[3]) +SMALL))

#define DIVBFLUXCTSTAGPRIM(p,i,j,k) ((                                  \
                                      DIVBDIFFFLUXCTSTAGPRIMx(p,i,j,k)/dx[1] + DIVBDIFFFLUXCTSTAGPRIMy(p,i,j,k)/dx[2] + DIVBDIFFFLUXCTSTAGPRIMz(p,i,j,k)/dx[3] \
                                      )*(DIVBNORMFLUXCTSTAGPRIM(p,i,j,k)))







////////////////////////////////////
///
/// FLUXCD DIVB
///
////////////////////////////////////

//////////////////////
/// CONSERVED quantity
/// FLUXCD divb lives at CENT
#define DIVBNORMFLUXCD(p,i,j,k) (MAX(MAX(dx[1],dx[2]),dx[3])/(MYGDET(i,j,k,CENT)*fabs( \
                                                                                      Fg(p,ip1mac(i),j,k,B1) + Fg(p,im1mac(i),j,k,B1) \
                                                                                      +Fg(p,i,jp1mac(j),k,B2) + Fg(p,i,jm1mac(j),k,B2) \
                                                                                      +Fg(p,i,j,kp1mac(k),B3) + Fg(p,i,j,km1mac(k),B3) \
                                                                                      )+SMALL))

#define DIVBFLUXCD(p,i,j,k)  (0.5*(                                     \
                                   (Fg(p,ip1mac(i),j,k,B1) - Fg(p,im1mac(i),j,k,B1))/dx[1] \
                                   +(Fg(p,i,jp1mac(j),k,B2) - Fg(p,i,jm1mac(j),k,B2))/dx[2] \
                                   +(Fg(p,i,j,kp1mac(k),B3) - Fg(p,i,j,km1mac(k),B3))/dx[3] \
                                   )*DIVBNORMFLUXCD(p,i,j,k))


//////////////////////
/// PRMITIVE quantity
/// FLUXCD divb lives at CENT
#define DIVBNORMFLUXCDPRIM(p,i,j,k) (MAX(MAX(dx[1],dx[2]),dx[3])/(MYGDET(i,j,k,CENT)*fabs( \
                                                                                          FgN(p,ip1mac(i),j,k,B1) + FgN(p,im1mac(i),j,k,B1) \
                                                                                          +FgN(p,i,jp1mac(j),k,B2) + FgN(p,i,jm1mac(j),k,B2) \
                                                                                          +FgN(p,i,j,kp1mac(k),B3) + FgN(p,i,j,km1mac(k),B3) \
                                                                                          )+SMALL))

#define DIVBFLUXCDPRIM(p,i,j,k)  (0.5*(                                 \
                                       (FgN(p,ip1mac(i),j,k,B1) - FgN(p,im1mac(i),j,k,B1))/dx[1] \
                                       +(FgN(p,i,jp1mac(j),k,B2) - FgN(p,i,jm1mac(j),k,B2))/dx[2] \
                                       +(FgN(p,i,j,kp1mac(k),B3) - FgN(p,i,j,km1mac(k),B3))/dx[3] \
                                       )*DIVBNORMFLUXCDPRIM(p,i,j,k))










////////////////////////////////////
///
/// SETTING MACROS for all divb methods
///
////////////////////////////////////


/// poles defined as divb=0, can't divide due to singularity (could use
/// volume regularization)
#define SETFDIVBFLUXCTTOTH(divb,U,i,j,k)     {DIVBCONDITION(U,i,j,k){ divb = fabs(DIVBFLUXCTTOTHPRIM(U,i,j,k)) ;} else divb = 0.;}
#define SETFDIVBFLUXCTTOTHPRIM(divb,p,i,j,k)     {DIVBCONDITION(p,i,j,k){ divb = fabs(DIVBFLUXCTTOTHPRIM(p,i,j,k)) ;} else divb = 0.;}


#define SETFDIVBFLUXCTSTAG(divb,U,i,j,k) {DIVBCONDITION(U,i,j,k){ divb = fabs(DIVBFLUXCTSTAG(U,i,j,k)) ;} else divb = 0.;}
#define SETFDIVBFLUXCTSTAGPRIM(divb,p,i,j,k) {DIVBCONDITION(p,i,j,k){ divb = fabs(DIVBFLUXCTSTAGPRIM(p,i,j,k)) ;} else divb = 0.;}


#define SETFDIVBFLUXCD(divb,U,i,j,k)     {DIVBCONDITION(U,i,j,k){ divb = fabs(DIVBFLUXCD(U,i,j,k)) ;} else divb = 0.;}
#define SETFDIVBFLUXCDPRIM(divb,p,i,j,k)     {DIVBCONDITION(p,i,j,k){ divb = fabs(DIVBFLUXCDPRIM(p,i,j,k)) ;} else divb = 0.;}

//#define SETFDIVB(divb,p,i,j,k) SETFDIVBFLUXCT(divb,p,i,j,k)


