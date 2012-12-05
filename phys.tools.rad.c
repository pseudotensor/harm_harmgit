#ifdef RADIATION
  ldouble eup[4][4],elo[4][4];
  pick_T(emuup,ix,iy,iz,eup);
  pick_T(emulo,ix,iy,iz,elo);
  ldouble Rij[4][4];
  calc_Rij(pp,Rij);
  boost22_ff2zamo(Rij,Rij,pp,gg,eup);
  trans22_zamo2lab(Rij,Rij,gg,elo);
  indices_2221(Rij,Rij,gg);

  //to move gdet in/out derivative:
  //here, up in metric source terms, in u2p and p2u, as well as in finite.c with del4[]
