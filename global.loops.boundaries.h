
/*! \file global.loops.boundaries.h
    \brief Boundary Conditions Loop definitions/macros

    // BOUNDARY CONDITION RELATED LOOPS
    
*/





//////////////////////////////////////
///
/// boundary condition loops:
///
/// The below loops have been generalized for any ORDERSTORAGE, which may help when boundary conditions are expensive for memory
//////////////////////////////////////

/// general loop for boundary conditions
#define BOUNDLOOPF(i,in,out) for(i=in;i<=out;i++)


/// for X1-dir [Note always uses innormal and outnormal so assumed to use this X1-dir loop first]
#define LOOPX1dir LOOPORDER1(,BOUNDLOOPF(j,innormalloop[2],outnormalloop[2]),BOUNDLOOPF(k,innormalloop[3],outnormalloop[3])) LOOPORDER2(,BOUNDLOOPF(j,innormalloop[2],outnormalloop[2]),BOUNDLOOPF(k,innormalloop[3],outnormalloop[3])) LOOPORDER3(,BOUNDLOOPF(j,innormalloop[2],outnormalloop[2]),BOUNDLOOPF(k,innormalloop[3],outnormalloop[3]))

/// for X2-dir [Note innormal and outnormal for 3rd direction and full inbound and outbound for 1st direction, so assume 1st came first and 3rd comes next]
#define LOOPX2dir LOOPORDER1(BOUNDLOOPF(i,inboundloop[1],outboundloop[1]),,BOUNDLOOPF(k,innormalloop[3],outnormalloop[3])) LOOPORDER2(BOUNDLOOPF(i,inboundloop[1],outboundloop[1]),,BOUNDLOOPF(k,innormalloop[3],outnormalloop[3])) LOOPORDER3(BOUNDLOOPF(i,inboundloop[1],outboundloop[1]),,BOUNDLOOPF(k,innormalloop[3],outnormalloop[3]))

/// for X3-dir [ Note only use of inbound and outbound, so assumes 1st and 2nd directions came first]
#define LOOPX3dir LOOPORDER1(BOUNDLOOPF(i,inboundloop[1],outboundloop[1]),BOUNDLOOPF(j,inboundloop[2],outboundloop[2]),) LOOPORDER2(BOUNDLOOPF(i,inboundloop[1],outboundloop[1]),BOUNDLOOPF(j,inboundloop[2],outboundloop[2]),) LOOPORDER3(BOUNDLOOPF(i,inboundloop[1],outboundloop[1]),BOUNDLOOPF(j,inboundloop[2],outboundloop[2]),)


/// loop over all directions (block of data in i,j,k only going out as far as the boundary conditions would go out)
#define LOOPXalldir LOOPORDER1(BOUNDLOOPF(i,inboundloop[1],outboundloop[1]),BOUNDLOOPF(j,inboundloop[2],outboundloop[2]),BOUNDLOOPF(k,inboundloop[3],outboundloop[3])) LOOPORDER2(BOUNDLOOPF(i,inboundloop[1],outboundloop[1]),BOUNDLOOPF(j,inboundloop[2],outboundloop[2]),BOUNDLOOP3(k,inboundloop[3],outboundloop[3])) LOOPORDER1(BOUNDLOOPF(i,inboundloop[1],outboundloop[1]),BOUNDLOOPF(j,inboundloop[2],outboundloop[2]),BOUNDLOOPF(k,inboundloop[3],outboundloop[3]))








//////////////////////////////////////
///
/// boundary loops for CENT quantities (these 1D loops need no generalization for any ORDERSTORAGE)
///
//////////////////////////////////////
///
/// normally start is a lower index than end and then loops take into account how to deal with this
///
/// most general loop ignorant of start begin before end
#define LOOPBOUNDGENMORE(i,iloopstart,iloopend,iloopstep) for(i=iloopstart;i<=iloopend;i+=iloopstep)
/// make start near boundary so periodic works for MAXBND up to 2N
#define LOOPBOUNDINGEN(i,start,end) for(i=end;i>=start;i--)
/// bound entire region inside non-evolved portion of grid
/// special horizon code still needed for multiple CPUs -- GODMARK -- unless chose range of CPUs inside horizon to apply inner BCs
//#define LOOPBOUNDINMOREGEN(i,start,end,ri) for(i=ri-1;i>=start;i--)
/// above "MORE" version is not more enough with "-1" and this technique is out-of-date compared to new grid sectioning approach
#define LOOPBOUNDINMOREGEN(i,start,end,ri) LOOPBOUNDINGEN(i,start,end)
#define LOOPBOUNDOUTGEN(i,start,end) for(i=start;i<=end;i++)

/// user should use set_boundloop() to control shifts
#define LOOPBOUND1INSPECIAL  LOOPBOUNDINMOREGEN(i,inoutlohi[POINTDOWN][POINTDOWN][1],inoutlohi[POINTDOWN][POINTUP][1],ri)
#define LOOPBOUND1IN  LOOPBOUNDINGEN (i,inoutlohi[POINTDOWN][POINTDOWN][1],inoutlohi[POINTDOWN][POINTUP][1])
#define LOOPBOUND1OUT LOOPBOUNDOUTGEN(i,inoutlohi[POINTUP][POINTDOWN][1],inoutlohi[POINTUP][POINTUP][1])

#define LOOPBOUND2IN  LOOPBOUNDINGEN (j,inoutlohi[POINTDOWN][POINTDOWN][2],inoutlohi[POINTDOWN][POINTUP][2])
#define LOOPBOUND2OUT LOOPBOUNDOUTGEN(j,inoutlohi[POINTUP][POINTDOWN][2],inoutlohi[POINTUP][POINTUP][2])

#define LOOPBOUND3IN  LOOPBOUNDINGEN (k,inoutlohi[POINTDOWN][POINTDOWN][3],inoutlohi[POINTDOWN][POINTUP][3])
#define LOOPBOUND3OUT LOOPBOUNDOUTGEN(k,inoutlohi[POINTUP][POINTDOWN][3],inoutlohi[POINTUP][POINTUP][3])

#define LOWERBOUND1 (inoutlohi[POINTDOWN][POINTDOWN][1])
#define UPPERBOUND1 (inoutlohi[POINTUP][POINTUP][1])



