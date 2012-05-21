
// whether metric mixes r-\phi or \theta-\phi.  Alows for optimizations since not often mixing with \phi
// Allow g_{\theta\phi} ?
#if(ALLOWMETRICROT==0)
#define DOMIXTHETAPHI 0 // choice, for optimizing
#else
#define DOMIXTHETAPHI 1 // NO choice, since rotation allows this mixing
#endif


