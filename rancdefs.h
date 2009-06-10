
// OPENMPMARK: below aa,S,n were static and calledranc was global, but ranc() may be called by multiple threads, so required to be global and use omp critical for that region.  Otherwise would have to make threadprivate and copyin() everytime think used, which is nasty.
int rancaa[NRANC];
int rancS[NRANC];
int rancvaln;
