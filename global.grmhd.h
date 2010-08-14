
#define DOINGGRMHDTYPECODE 1 // always 1
#define DOINGGRRAYTYPECODE 0 // always 0
#define DOINGLIAISONTYPECODE 0 // always 0
#include "global.general.h"

/*
  Generally, switches that are performance related are here, while others are in init with variable complements for flexibility during runtime.

  Not everything is in here so can compile quicker on simple changes.

*/




/////////////////////////////////////////
//
// Comments on comments
//
//////////////////////////////////////////

//
// CHANGINGMARK: Most recent changes or issues
//
// GODMARK: Generally needs consideration but may not be crucial
//
// SUPERGODMARK: Needs urgent attention or workaround
//
// SUPERMARK: Something changing
//
// SASMARK: Changes by Sasha
//
// SECTIONMARK: Issue/Comment with regards to grid sectioning
//
// TESTMARK: Something was/is for testing purposes only
//
// NEWMARK: New code in development
//
// **TODO** : Things to do in code
//
// GLOBALMARK: Location where global positional array is used
//
// ISSUESMARK: List of known issues of importance
//
// USAGEMARK: Notes on usage of code
//
// OPTMARK: Comment about performance optimizations
//
// SUPERNOTE: Important note should consider anytime in that file
//
// OPENMPNOTE: Note related to OpenMP
//
// OPENMPOPTMARK: Optimization note or suggestion
//
/////////////////////////////////////////



////////////////////////////////////
//
// Debugging comments
//
////////////////////////////////////

//
// debugging on Harvard's Bluegene:
// Compile with -g -O0
// On BG/P, core files are text files.  Look at the core file with a text editor, focus on the function call chain; feed the hex addresses to addr2line.
// addr2line -e your.x  hex_address
// 
// tail -n 10 core.511 | addr2line -e your.x
// 
// Use grep and word-count (wc) to examine core files :
//  grep hex_address “core.*” | wc -l
// 
// You can get the instruction that failed by using objdump:
// powerpc-bgp-linux-objdump -d your.x >your.dump
// You can locate the instruction address in the dump file, and can at least find the routine where failure occurred, even without –g.
// 
// If your application exits without leaving a core file,
// set the env variable BG_COREDUMPONEXIT=1

/* coreprocessor.pl … coreprocessor perl script in : */

/* /bgsys/drivers/ppcfloor/tools/coreprocessor */

/* online help via :  coreprocessor.pl –help */

/* Can analyze and sort text core files, and can attach to hung processes for deadlock determination.   */

/* Normal use is “gui” mode, so set DISPLAY then: */

/* coreprocessor.pl –c=/bgusr/username/rundir –b=your.x */

/* click on “Select Grouping mode”,  */
/* select  “Stack Trace (condensed)” */
/* click on source statement of interest */

/* There is also a non-gui mode: coreprocessor.pl –help */


// Ensure to test for external declaration consistency using checkexterns.sh.  See end of global.general.h





////////////////////////////////////
//
// Performance Notes
//
////////////////////////////////////

// 0) Consider file writing and other bottlenecks
//
// PRODUCTION==0 vs. 1 changes performance by about 20% with debugfail==2 set
// E.g.: Type of problem can change performance.  For example, same 64x32 model with R0=0 Rout=40 gets 70K ZCPS while R0=-8.3 and Rout=1E10 gets 62K ZCPS.  Both have grid sectioning.  The change in performance here occurs because one has more failures than the other.  Set PRODUCTION 0 allows no file writing of those failures and gets that 62K ZCPS run back to 70K.
//
// DODIAGS==0 or 1 can change things alot due to per-substep diagnostics being enabled.
// E.g.: File writing in general should be stretched in period to avoid slowing down code.  File writing can severely slow things down if too frequent.  Consider ROMIO and Jon's non-blocking file writing.
// E.g.: Same model described above w.r.t. PRODUCTION=0/1 goes to 80K ZCPS with DODIAGS=0.  Ensure that ENERDUMPTYPE (and image/etc. creation) does not have too small a period.
//
//
//
//
// 1) Consider cache use.
//
// For example, on 4 core system (Intel Core2), found that 512x32 model runs at 62K ZCPS while otherwise similar 64x32 model runs at 64K ZCPS.  As below states, could fit about 100 zones of data into cache, and having N1=64 (or N1M=72) allows entire line of data to fit into cache.  This reduces cache misses.
//
// On the other hand, if doing grid sectioning with large number of radial zones (e.g. 1024), then optimal to have entire N1M line on one node.  This will necessarily cause additional cache misses, but otherwise nodes won't be used and the performance drop is factors of 2X or more instead of 5% as above.
//
//
// 2) Consider each system's cache.
//
// e.g. TACC Ranger uses 4 sockets of Barcelona CPU:
// http://images.anandtech.com/reviews/cpu/amd/phenom2/barcelona-block-diagram.jpg
//  Each core has 64KB private L1 cache (32K L1I and 32K L1D) (I=insruction D=data)
//  Each core has *private* 512KB L2 cache.  This is nice compared to the shared L2 cache of my ki-rh42, even if that cache is 2X larger.
//  Each Socket of 4 cores has *shared* 2MB L3 cache.
//  Each node has 4 sockets each with 4 cores = 16 cores/node
//  Node has 16 cores with 32GB/node
//  The memory bus runs at 533Mhz with 2 channels total for *all* 16 cores!
// So staying in cache is *critical* since otherwise all 16 cores compete for 2 channels of memory.
// If using OpenMP, then one avoids otherwise MPI-excessive boundary cells, so good. [So this really makes OpenMP save on memory!!! Very important!]
// Can't fit entire problem in cache.  In reality, need to know how much memory use over longer periods of time.  Perhaps think per-line.  Then only need to ensure that can do a single line.  Then can do up to about N1=160 total on those 16 cores before cache-misses will occur.
//
// TACC Lonestar has 4MB L2 cache with 2.66Ghz Intel Xeon 5150.  TACC page calls it "smart" L2 cache.
// Appears to have private 2MB per core.
// http://www-dr.cps.intel.com/products/processor/xeon5000/specifications.htm?iid=products_xeon5000+tab_specs
// http://services.tacc.utexas.edu/index.php/lonestar-user-guide
// http://processorfinder.intel.com/List.aspx?ProcFam=528&sSpec=&OrdCode=
// http://processorfinder.intel.com/details.aspx?sSpec=SL9RU
// http://www.linuxdevices.com/files/misc/intel_5100.jpg
// http://www.hardwarezone.com/img/data/articles/2006/2002/Bensley-block-diagram-2.jpg
//
// Lonestar + Ranger -> Ranch for data storage:
// http://services.tacc.utexas.edu/index.php/ranch-user-guide/
// tar cvf - thickdisk1 | ssh ${ARCHIVER} "cat > ${ARCHIVE}/thickdisk1.tar"
//
// With ssh (according to TACC):
// cat myfile.tar | ssh target_machine "cd /target_dir/; tar xvf - "
//or, for transfers between TACC systems, you could use the version of
//gsissh that we have installed, which allows you to do the handshake
//encrypted but send the data unencrypted:
// cat myfile.tar | gsissh -oNoneEnabled=yes -oNoneSwitch=yes target_machine "cd target_dir; tar xvf - "
//The overall consensus here is that the fastest way to do what you want
//is to expand the files to $SCRATCH in Ranger and then bbcp them down
//to your machine. This should be faster than using the ssh/cat/tar that
//you are using now.
//
// other options:
// http://moo.nac.uci.edu/~hjm/HOWTO_move_data.html
//
// bbcp much faster:
// http://www.slac.stanford.edu/~abh/bbcp/
// http://pcbunn.cithep.caltech.edu/bbcp/using_bbcp.htm
// bbcp -P 5 -k -a ....
//
// On jdexter@@luigi : getting files from ki-rh42:
// ~/bin/bbcp -r -P 5 -a -k -T 'ssh -x -a %I -l %U %H bbcp' jon@ki-rh42.slac.stanford.edu:/media/disk/jon/thickdisk1/dumps/fieldline* jon@ki-rh42.slac.stanford.edu:/media/disk/jon/thickdisk1/dumps/dump0000* .
//
// On QB to ranch (problem with -a into ranch -- have to manually  make directories used for recovering -- stupid):
// ~/bin/bbcp -r -P 5 thickdisk5fghij tg802609@ranch.tacc.utexas.edu:
//
// If copy failed for whatever reason (QB went down), then have to copy per directory since -a -k don't work properly with -r.  So one has to complete:
//~/bin/bbcp -a -k -P 5 * tg802609@ranch.tacc.utexas.edu:thickdisk5fghij/
// cd dumps/
//~/bin/bbcp -a -k -P 5 * tg802609@ranch.tacc.utexas.edu:thickdisk5fghij/dumps/
// cd ../images/
//~/bin/bbcp -a -k -P 5 * tg802609@ranch.tacc.utexas.edu:thickdisk5fghij/images/
// If get message of append couldn't be completed, delete that file on ranch
// If get message that not enough space to copy files, then contact TACC.
//
// bbcp -a will not copy if prior copy failed in some way that bbcp thinks the file is changed.  The file will be skipped!
// Can create a file list:
//
// ls > filelist.txt
// OR to ignore some files:
// ls --ignore=dumps --ignore=images --ignore=thickdiskr3 --ignore=filelist.txt --ignore=err.txt > filelist.txt
//
//
// and copy and check err.txt file for such messages "append not poosible":
// ~/bin/bbcp -a -k -P 5 -I filelist.txt -l err.txt tg802609@ranch.tacc.utexas.edu:thickdiskr2/dumps/
//
// OR to avoid time-out issues with large file lists:
// ensure -b value is at least twice -s value
// can make s up to (say) 32 or 64 on LAN, but keep to order 4-10 on WAN
//
//////// ~/bin/bbcp -s 4 -b 8 -w 2m -a -k -P 5 -I filelist.txt -l err.txt -T 'ssh -x -a -oConnectTimeout=0 -oFallBackToRsh=no %I -l %U %H bbcp' -S 'ssh -x -a -oConnectTimeout=0 -oFallBackToRsh=no %I -l %U %H bbcp' tg802609@ranch.tacc.utexas.edu:thickdiskr2/dumps/
//
// ~/bin/bbcp -s 4 -b 8 -w 2m -a -k -P 5 -I filelist.txt -T 'ssh -x -a -oConnectTimeout=0 -oFallBackToRsh=no %I -l %U %H bbcp' -S 'ssh -x -a -oConnectTimeout=0 -oFallBackToRsh=no %I -l %U %H bbcp' tg802609@ranch.tacc.utexas.edu:thickdiskr2/dumps/
//
// OLD: Best not use bbcp's -r for recursive since won't recover properly after using -a -k when copy failed.
//
//
//
//
// Edited the bbcp code in bbcp_FileSpec.C to overwrite file if append not (sic) poosible, but on ranch can't compile.
//
// Note that if use -l err.txt, then may not tell you if completed successfully or not!  Just stops.  While NOT using -l err.txt means don't know if had any "poosible" problems.
//
// In the end, this worked:
// 1) use bbcp as below.  Note that bbcp is the only program I've found that resumes.  For example, for one directory do:
//
// ls --ignore=err.txt --ignore=filelist.txt grmhd > filelist.txt
// [grmhd included above because otherwise might have * in name and bbcp doesn't like that]
// -b + 100 helps avoid stalls with ranch and tape processing on ranch side.
//
// ~/bin/bbcp -a -k -f -b +100 -P 5 -V -I filelist.txt -T 'ssh -x -a -oConnectTimeout=0 -oFallBackToRsh=no %I -l %U %H bbcp' -S 'ssh -x -a -oConnectTimeout=0 -oFallBackToRsh=no %I -l %U %H bbcp' tg802609@ranch.tacc.utexas.edu:thickdiskr3.abcdefgh/dumps/ &> err.txt
//
// Can also try without -I filelist.txt and use recursive copy in new version:
// Add: -r  and give path name(s) just before destination
// E.g.:
//
// ~/bin/bbcp -b +100 -a -k -f -r -P 5 -V -T 'ssh -x -a -oConnectTimeout=0 -oFallBackToRsh=no %I -l %U %H bbcp' -S 'ssh -x -a -oConnectTimeout=0 -oFallBackToRsh=no %I -l %U %H bbcp' thickdisk14.ab thickdiskrr2.ab tg802609@ranch.tacc.utexas.edu: &> err.txt
//
//
// 2) Copy each directory (the root data, dumps, images).  Can be done at the same time in different shells.  Don't use recursive copy in bbcp since won't resume correctly.
//
// 3) Once copy "done", need to check if really copied and if all files copied and if correctly copied.
//
// 4) Repeat copy and ensure stderr (stored in err.txt) never outputs "poosible" or stops copying certain files.  This checks that append and all things are making sense to bbcp.  bbcp can append incorrectly and make files too short and even too long!
//
// 5) Repeat copy and output stderr to file so can:
//    grep "already been copied" err.txt | wc
//    This checks number of actual files copied.  bbcp can skip files in the list provided, or even if looping on shell, or if using shell expansion.  So unreliable.
//
// 6) Check result on remote machine per directory:
//
// ls --ignore=dumps --ignore=images --ignore=filelist.txt --ignore=err.txt |wc
//
// 6.5) Check filelist.txt to ensure actually copied all files:
//
// wc filelist.txt
// OR for recursive folders:
// find <list of space separated folders> > filelist.txt | wc
// E.g.:
// find thickdiskrr2.ab thickdisk14.ab > postfilelist.txt
//
//
// 7) If numbers match, then likely have copied everything properly.  Should be correct size too, but can check manually that files have uniform looking sizes.  du -s doesn't help much since on ranger sizes are measured oddly.
// 
//
// Tried uberftp and globus-url-copy, but problems
//
//uberftp -binary -cksum on -keepalive 60 -parallel 10 -resume thickdiskr2 ranch.tacc.utexas.edu
//then enter and do:
// resume thickdiskr2
// put -r thickdiskr2
// for recursive copy.  Not sure how failure resume works if at all.
//
// Problem is that certificate expires very soon and kills copy, so useless to copy large files.
//
// Also tried:
// globus-url-copy -b -r -restart -p 5 file://thickdiskr2  gsiftp://ranch.tacc.utexas.edu:
//
// But failed to figure out how to get it's special certificate using that other program with a similar name.
//
//
//
//
// Teragrid LONI QueenBee (Queen Bee):
// http://www.loni.org/systems/
// http://www.loni.org/teragrid/users_guide.php
// 1) Connect to (e.g.) abe: ssh jmckinne@abe.ncsa.uiuc.edu
// 2) Setup proxy: myproxy-logon -l jmckinne  -s myproxy.teragrid.org
// 3) Use NCSA Portal Password
// 4) Connect to QB: gsissh login1-qb.loni-lsu.teragrid.org

// QUEENBEE: Must set GETTIMEOFDAYPROBLEM 1

//
//
// Teragrid NCSA Abe:
// http://www.teragrid.org/userinfo/hardware/resources.php?select=single&id=50&PHPSESSID=0
// http://www.dell.com/content/products/productdetails.aspx/pedge_1955?c=us&l=en&s=corp
// 5300 quad-core Xeon 5000 series
// 2*4MB L2 cache (shared like ki-rh42).
//
// tgusage -u jmckinne
// tgusage -a TG-AST080026N
// tgusage -a TG-AST080025N
//
//
// Intel Core2 has same(?) as Barcelona but no L3 cache.
// JCM's system is: Intel(R) Core(TM)2 Extreme CPU Q6800  @ 2.93GHz stepping 0b
// Appears to have 2x4MB=8MB (shared) L2 cache and 32K/32K D/I L1 caches
// This 8MB of L2 cache is (however) shared so that 2 cores share 4MB and the other 2 shared 4MB.  So unless tune which core gets which process, each use of a new core slows down the L2 cache bandwidth by roughly 2X.  May be consistent with performance results below.
// So 2X cache compared to Ranger.
// http://www.intel.com/design/processor/datashts/316852.htm
//
// Intel Core i7:
// http://en.wikipedia.org/wiki/Intel_Core_3
// Note it has 32KB L1 D cache per core
// 256KB L2 cache (I+D) per core
// 8MB L3 cache (I+D) shared between *all* cores.
// Tri-channel memory
//
//
//
///////////////////////////////////////////////////////
// Testing memory contention and whether sit entirely within L2 cache on JCM's system without dual-channel access currently (so strong test of whether run fits in L2 cache or not):
//
// Multiple runs (performance for each core):
//               #Cores:        1     2     3       4
// MAXBND==2 DONOR    64x32:   60K                 27K-30K
// MAXBND==4 PARALINE 64x32:   50K   45K           24K-27K
// MAXNBD==2 DONOR    4x2:     30K   30K   30K     30K       [so finally fit entirely within L2 cache]
// MAXNBD==4 PARALINE 4x2:     18K   18K   18K     18K       [so finally fit entirely within L2 cache]
// MAXBND==4 PARALINE 4x8:     30K   29.5K  29.3K  29K       [so finally fit entirely within L2 cache]
// MAXBND==4 PARALINE 16x8:    50K   48K    48K    48K       [so kinda hits memory]
// MAXBND==4 PARALINE 32x16:   56K   55K    36K    31K       [so 32x16*(2cores) is ok in cache, but not more (i.e. not this run with 32*16*3 or 32*16*4)]
// This suggests an optimal OpenMP choice is N1*N2*N3=32*32 per *node* to avoid memory contention and stay within L2 cache.  For elongated cells this is 512x2 in 2D or 256x2x2 in 3D.  But small cells leads to inefficiency in extra computations near boundary cells.  So probably want 256x4 in 2D and 64x4x4 in 3D.  On Ranger have to take (1/2) of this since L2+L3 is half the size!  So on Ranger use 128x4 in 2D and 32x4x4 in 3D with fit into L2+L3 cache.
//
//
// For above, probably good 4 core performance because staying in *L1* cache mostly, not L2.  This is consistent with below results.
// Then hit is because ki-rh42's 4 cores feed off 2 L2 cache's so that bandwidth is half when 3-4 cores are running.
//
///////////////////////////////////////////////////////
//
//
///////////////////////////////////////////////////////
// New code with even more optimizations (esp. fluxct and no eomfunc[NPR] data and now also reduction of symmetric matrices):
//
// Multiple runs (performance for each core):
//               #Cores:        1     2     3       4
// MAXBND==4 PARALINE 16x8:    51K   51K    51K    50K       [4 core perf not helped with symm. matrix reductions]
// MAXBND==4 PARALINE 32x16:   60K   60K    54K    50K       [so cache misses not as bad with new code. Reducing symmetric matrices increased performance at 4 cores by 25%, so good.]
// MAXBND==4 PARALINE 64x16:   60K   59K    48K    36K       [sym matrix fix increased perf for 4 cores by 25%]
//
//
// OpenMP (i.e. USEOPENMP==1 in makehead.inc) tests:
//
//               #Cores:        1     2     3       4
// MAXBND==4 PARALINE 16x8:    37K   48K   58K     66K       [horrendous 4 core performance -- 10% improvement with minchunk of 10 or static schedule compared to guided.  Didn't help with 1 core overhead.  Reducing memory overhead didn't help, suggesting it's all OpenMP overhead!]
// MAXBND==4 PARALINE 64x16:   48K   78K   88K     99K       [Unsure what's memory vs. OpenMP overhead, but big hit compared to above runs!  Reducing memory overhead didn't help at all for 4 cores!  Suggests all OpenMP overhead!  Removed as much of private() [all essentialy] and copyin() [as much as I could] and little change! Even ALLOWKAZEOS is off! -- Must be EOS ptr functions?  Just stupid addresses!]
// MAXBND==4 PARALINE 32x16:   45K   69K   83K     92K
//
// Even with memory contention, 4 cores here should get about 4*40K=160K or 4*30K=120K.  Appears to be 20%-25% overhead per core no matter how many (even 1) cores.
//
// Also appears to be extreme overhead no matter what when doing <200 or so iterations.  Requires per node memory be (say) 64x16 or larger.  But L2 cache requires 64*16*2 or smaller.  So seems 64x16 - 64x32 (or 32x32) is sweet spot for minimizing OpenMP overhead and L2 cache misses.
//
//
// To test OpenMP overhead, I commented out all #pragma's but left -openmp during compilation and checked 1 core:
//
//               #Cores:        1     2     3       4
// MAXBND==4 PARALINE 64x16:   61K
//
// So not -openmp itself, but #pragma's cause the problem.
//
//
// Now I tried only commenting out inner loop pragmas:
//
//               #Cores:        1     2     3       4
// MAXBND==4 PARALINE 64x16:   46K
//
// So not -openmp or #pragma omp for, but #pragma omp parallel is problem.
//
//
// Now I tried commenting out loops and choosing OPENMPGLOBALPRIVATE (etc.) to be nothing (can do this since globals preserved on master):
// Also commented out #pragma omp threadprivate....
// And commented out parallel copyin( -> //copyin(
//
//               #Cores:        1     2     3       4
// MAXBND==4 PARALINE 64x16:   61K
//
// So problem appears to be copyin() .... with same 64x16 test:
//
// TEST1CORE: without copyin(EOS ptr's) reach 52K for 1 core, still overhead?
// TEST1CORE: without copyin(EOS ptrs' + all loop stuff except ploop): reach 52K : still overhead!
// TEST1CORE: without copyin(EOS ptrs' + all loop stuff): reach 57K : finally little overhead!
//
// NEWCODE1: with all required things as thread private and copyin() [had too many things for loops]: 51K
// NEWCODE2: with EOS ptrs totally fixed (outside parallel and not thread private anymore): 54K-55K for 1 core [good!] 102K or up to 112K for 4 cores [still bad]
// NEWCODE3: with no {ijk}curr, whichuconcalc, etc. that required extensive code adjustments for the EOS-related stuff to force modularity: 57K for 1 core and 104K for 4 cores : So quite little OpenMP overhead now for 1 core, but still there for 4 cores.
//
// So 1 core overhead is quite small.  Now unsure why 4 core bad since not memory limited since changing cores doesn't change problem size.
//
//
// Installed all memory into ki-rh42 (8GB total) so uses dual-channel memory access:
// 1) OpenMP 4 cores went from 86K -> 107K for 512x32 testreality problem.
// 2) MPI 4 cores (no grid sectioning so uses all cores) went from 90K -> 110K
// So despite all Jon's hard work with OpenMP, MPI still faster?  Can't be large overhead with OpenMP at 512x32.
// Overall, performance is nearly as expected for simulation that doesn't fit entirely into L2 cache that would be roughly 31K/core * 4 = 120K. 
//
//
//  
//                   #Cores:        1     2     3       4
// MPI MAXBND==4 PARALINE 64x16:                       72K (64x16 total with 32x8 per core of 4 cores) -- only runs at 40% of core per process!
// OMP MAXBND==4 PARALINE 64x16:                       130K (only a bit over 2X expected performance from 57K/core from NEWCODE3.  Not sure what's limiting since should be all in L2 cache!)
//
// Clearly for small problem sizes per core, OpenMP is much more efficient than MPI.  For large problem sizes per core (e.g. 256x16 per core for 512x32 total), MPI is just slightly more efficient.  Thus, in general one is benefited by using OpenMP.  Unclear how things scale with many more processores or cores.
//
// Still problem is that efficiency from 1 core to 4 cores is 50% for either case, and again not because of L2 cache [Well, ki-rh42's L2 cache is effectively shared, so may explain it.]
//
// Suggests that, for Ranger at least, should try to fit into private 512KB L2 cache.  Can't use above tests to conclude what size this is for ki-rh42.  Can only run 2 cores and see when adding 2nd core slows down (assumes each core uses essentially its own nearby cache).  Once 2nd core slows things down significantly, then must have gone to memory.
//
// 2CORE ki-rh42 test (latest NEWCOD3) for fitting into L2 cache:
// Multiple runs (performance for each core, OpenMP disabled):
//               #Cores:        1     2    3     4
// MAXBND==4 PARALINE  64x16:   63K   63K
// MAXBND==4 PARALINE  64x32:   61K   60K
// MAXBND==4 PARALINE  64x64:   56K   58K
// MAXBND==4 PARALINE 128x128:  61K   57K             [Starting to hit memory]
// MAXBND==4 PARALINE 256x256:  60K   56K 50K    40K  [Starting to hit memory]
//
//
// OpenMP:
//               #Cores:        1     2    3     4
// MAXBND==4 PARALINE 256x256: 59K   95K  114K  116K  [Consistent with (roughly) 2 cores sharing 4MB of L2 cache, so with 3-4 cores goes much slower per core due to cache contention.]
// MAXBND==4 PARALINE  64x16:  59K   99K  122K  128K  ["" -- roughly]
// MAXBND==4 PARALINE  32x16:  59K   90K  115K  119K  [Linux probably doesn't schedule processes based upon L2 cache association with cores, so can flip around and almost as bad as all 4 cores sharing L2 cache.]
//
//
// Perf with STAG+DISS+DISSVSR+LUMVSR    = 30K -- so about 2X slower.
// Perf with STAG+DISS+DISSVSR+LUMVSR+3D = 11K -- so about 6X slower.
//
// TODO:
// 1) if limiting interpolation (e.g. for stag or rescale in stag), then pass that fact rather than using global npr2interp.  Ensure all interior loops use that passed data rather than globals.
// 2) For KAZ stuff, probably fine.
//
//
// Performance on Ranger (All MAXBND==4 PARALINE) for 1000 steps with DODIAGS=0 and PRODUCTION 1 and TIMEORDER=2 and FLUXB=FLUXCTTOTH and DODISS=DODISSVSR=DOLUMVSR=0:
//
// N1xN2    #NODES  #OPENMP/task  #MPItasks  ncpux?      PERF   Eff
//
// 64x64:      1       1             1        1x1x1       35K  100%
// 64x64:      1       4             4        2x2x1      403K   72%
// 64x64:      1       1            16        4x4x1      407K   73%
// 64x64:      1      16             1        1x1x1       71K   13%
// 64x64:      1       8             2        2x1x1      215K   39%
// 64x64:      2       1            32        8x4x1      940K   84%
// 64x64:      2       4             8        4x2x1      950K   85%
// 64x64:      2       1            64        8x8x1     1793K   80%
// 64x64:      2       4            16        4x4x1     1862K   83%
// 64x64:      2       1           256        16x16x1   7098K   80%   
// 64x64:      2       4           64         8x8x1     6766K   76%

//
// 64x16:      1       1             1        1x1x1       35K
// 64x16:      1       4             4        2x2x1      406K   72%
// 64x16:      1       1            16        4x4x1      471K   84%
// 64x16:      1      16             1        1x1x1       55K   10%
// 64x16:      1       8             2        2x1x1      191K   34%
// 64x16:      2       1            32        8x4x1      907K   81%
// 64x16:      2       4             8        4x2x1      777K   70% [repeated run got same result!?]
// 64x16:      2       1            64        8x8x1     1730K   77%
// 64x16:      2       4            16        4x4x1     1488K   66%
//
// So each Ranger core is almost 2X slower than ki-rh42.  That's AMD vs. Intel for you!
// So clearly bad to cross on PCI bus with memory as OpenMP has to when more than 4 threads with 1 thread per core.
// MPI seems to be doing fine at 64^2 on one node.
// Unclear why OpenMP is actually slower even for 64x16, which worked better on ki-rh42 by 2X!
//
//
//
// Performance on Ranger (All MAXBND==4 PARALINE) for 1000 steps with DODIAGS=0 and PRODUCTION 1 and TIMEORDER=2 and FLUXB=FLUXCTSTAG and DODISS=DODISSVSR=DOLUMVSR=1:
// Default: module unload mvapich2 pgi ; module load mvapich2 intel mkl
//
// N1xN2    #NODES  #OPENMP/task  #MPItasks  ncpux?      PERF   Eff
//
// 64x32x8:    1       1             1        1x1x1        6K
// 64x32x32:   1       1             1        1x1x1        8.7K        [to be used as reference for efficiency for 64x32x32 per MPI task runs]
// 64x32x8:    4       4             4        2x2x1       41K   43%
// 64x32x8:   16       1            16        4x4x1       84K   88%
// 64x32x8:   64       4            64        8x8x1      617K   40%
// 64x32x8:   64       4            64        8x8x1      624K   41%  [used purely static schedule, no user chunking]
// 128x64x8:  64       4            64        8x8x1      713K   46%  [used purely static schedule, no user chunking]
// 64x32x32:  64       4            64        8x8x1      884K   40%-57%  [used purely static schedule, no user chunking] [smaller % is when using correct reference point]
// 64x32x32:  64       4            64        8x8x1      896K   40%-58%  [Changed to mvapich/1.0.1 during compile and in batch script]
// 64x32x32:  64       4            64        8x8x1      888K   39%-58%  [Changed to openmpi/1.3 during compile and in batch script]
// 64x32x8:  256       1           256        8x8x4     1278K   83%
//
//
// On ki-rh42 with new code:
// N1xN2    #NODES  #OPENMP/task  #MPItasks  ncpux?      PERF   Eff
// 64x32x8:    1       1             1        1x1x1       17K
//
// From 11K -> 17K -- not that big a difference after a day of pain.
// 
//
// On Ranger: 4x4x4 per CPU with PARALINE+STAG+DISS+LUM (but R0=0 and Rout=10) on 2048 processors and 16x16x8 for ncpux?: late-time with 54K ZCPS with average of 10% fractional diagnostics
//
//
//FULLSTAG=UNFUDDLE,STAG, RK2, DODISS, currents, entropy evolution, etc.etc.
//HALFSTAG=UNFUDDLE,STAG,RK2, no features
//HALFTOTH=UNFUDDLE,Toth,RK2, no features
//ORIG=Original very lean (and crashy-unstable) 2D HARM, which is on unfuddle as "origcode"
//
//system   tile_size   cores   zcps   efff         CODE-MODE
//----------------------------------------------------------------------
//3D:
//lonestar  34x32x8    1        13K       1            FULLSTAG
//lonestar  34x32x8    1        14K       1            HALFSTAG
//lonestar  34x32x8    1        28K       1            HALFTOTH
//lonestar  34x32x8    1024     9.5M      71%          FULLSTAG
//lonestar  34x8x8     1        ?         1            FULLSTAG
//lonestar  34x8x8     512      3.3M      >50%?        FULLSTAG
//
//ki-rh42    32x16x8    1         14K      1           FULLSTAG
//ki-rh42    32x16x8    1         15K      1           HALFSTAG
//ki-rh42    32x16x8    1         28K      1           HALFTOTH
//ki-rh42    16x32x32   1         30K      1           HALFTOTH (Noble setup, who got 36K)
//2D:
//ki-rh42    64x64      1         57K       1          HALFTOTH
//ki-rh42    256x256    1         61K       1          HALFTOTH
//ki-rh42    64x64      1         80K       1          ORIG
//ki-rh42    256x256    1         84K       1          ORIG
//
//
//
///////////////////////////////////////////////////////
//
//
// http://www.intel.com/design/core2XE/documentation.htm
//
// 1) Best if arrays are powers of 2 so compiler bit shifts to index array
// 2) Try to keep all memory local.  Cache-lines are grabbed along fastest memory portions, so avoid loops that skip over indices.  Also best to compact multiple arrays together into a single per-point array or structure if those arrays are accessed together since then avoids cache miss when having to grab other array that will be displaced due to its full size.
// 3) Try to reduce memory footprint (total arrays, global or not) (and also reduce N1,N2,N3) so code+data fit into L2 cache (~4MB+)
// 4) Try to ensure code+data at function level can fit into L1 cache (32K)
// 5) Esp. for multi-core apps, want to fit MUCH of data into L2 cache
//    This corresponds to (currently) using 25K/cell including BZones
//    So (e.g.) N1M=N2M=12 (i.e. N1=N2=4) would fit completely into L2 cache of 4MB
//    In reality, just need to fit code+data along several cache lines to avoid memory being accessed.  Code is currently directed in memory along N3, then N2, then N1.  So optimal for N3>>N2>>N1.
//    Should allow N1=phi, N2=theta, N3=r in case want many r cells per core
//    Instead of changing [N1M][N2M][N3M], should change from [i][j][k] -> [k][j][i] for example.  So macrofy this.
//    Similar, but more involved than what Sasha does with his init's.
//
//
//
// http://www.eventhelix.com/RealtimeMantra/Basics/OptimizingCAndCPPCode.htm
//
// When arrays of structures are involved, the compiler performs a multiply by the structure size to perform the array indexing. If the structure size is a power of 2, an expensive multiply operation will be replaced by an inexpensive shift operation. Thus keeping structure sizes aligned to a power of 2 will improve performance in array indexing.
//
//
//  3) Note that for very small N1,N2,N3, ZCPS will appear to drop even if internal performance is the same.  This is because (e.g.) if N1=32 and N2=2, then in 2-direction it's really like doing N2=4 because of the extra 2 zones cmputed for the surface flux.  In reality not all parts of the code do this, so affected performance is between 1-2X factor.
//
// 

// Use pfmon (perfmon2) to look at cache and memory problems.
// Note that perfmon2 by default only looks at a single thread, but shows OpenMP usage for that thread so still useful.  Also shows pow() and other library calls that gprof does not.
// Note that perfmon2's HALT cycles may be larger for faster code in a repeated way (not just other procs eating time).  This doesn't mean if time without perfmon that will be slower.  For example, I found OpenMP'ing poledeath() led to more HALT cycles but clearly shorter wall time.
//
// See installperfstuff.txt and compilekernel.txt

// Use gprof:
// 1) Compile with -g -pg
// 2) Run test code
// 3) Run: gprof ./grmhd > prof.txt
// 4) Run: gprof -l ./grmhd > profbyline.txt
// 5) Look at prof.txt top and bottom portions to get idea of what is bottleneck
// 6) Use profbyline.txt to see if issue with bottleneck happens to be single lines that can be improved.

// In any case of pfmon or gprof, compile with -fno-inline (remove other inline commands) to see inlined functions expanded so can tell what interior functions take time within a function that was previously inlined.  Use __inline keyword in front of those functions one wants to see not inlined.

// Note that while gprof seems to miss some functions (e.g. pow() in icc library), pfmon catches it as pow.L .  

// One can use gprof-helper.c to have gprof measure multi-thread performance to see overhead.
// See inside of gprof-helper.c:
// do:
// 1) gcc -shared -fPIC gprof-helper.c -o gprof-helper.so -lpthread -ldl
// 2) (e.g.) LD_PRELOAD=./gprof-helper.so grmhd 4 1 1 1

// Note that both gprof and pfmon or any performance timer will slow down code.  Remove -pg from makefile to avoid slowdown as well.

//
//
//If not using the SIMD (single instruction, multiple data), or SSE type operations (which automatically are tried if using most compilers with optimizations), then roughly:
//  CPU
//  cycles | operation or procedure
//---------------------------------------------
//1          addition, subtraction, comparison
//2          fabs
//3          abs
//4          multiplication
//10         division, modulus
//20         sqrt
//50         exp
//60         sin, cos, tan
//80         asin, acos, atan
//100        pow
//10         Miss L1 cache
//50         Miss L2 cache if L3 cache present
//150        Branch misprediction
//200        Miss L3 cache or miss L2 cache if no L3 cache present
//200-1000+  Page fault
//
//Also, note that there are more than just math operations to worry about.  Another rule of thumb is that:
//
//page fault >>>>>> cache miss >> branch mistaken >> dependency chain >> non-sse-sqrt, division, and modulus > sse-sqrt etc. > everything else
//
//with ">>" being much greater than (as in 3-20 times).
//
//
////////////////////
//
// To get size of code's array elements:
//
// Use program "nm" to list objects and symbols, or use objdump -axfhp
// Use "nm -S --size-sort" to list by memory size

// To check global symbol sizes use:
//      nm -S --size-sort grmhd 
// or:  objdump -axfhp grmhd
//
// e.g. to check total size used do:
// 1) nm -S --size-sort grmhd | awk '{print $2}' > list.grmhd.sizes.txt
// 2) total=0
// 3) for fil in `cat list.grmhd.sizes.txt` ; do total=$(($total+0x$fil)) ; done
// 4) echo $total
// Now take that number and divide by N1M*N2M*N3M and that's approximate per zone bytes used  (must include boundary cells)
//
// Seems to be about 1250 elements/zone



// Checking FLOPS with papi
// http://www.cisl.ucar.edu/css/staff/rory/papi/flops.html
//
// 
// with perfmon2 (pfmon):
// pfmon --show-time -e FP_COMP_OPS_EXE ./grmhd.makedir.testpfmon_simple 1 1 1
// Note that X87_OPS_RETIRED:any doesn't seem to make sense
// Then take resulting number and divide by the number of seconds for running the code WITHOUT pfmon!  That's FLOPS.
// I get about 580MFLOPS  per core \sim 0.6GFLOPS per core
// Theoretical maximum is about 5.4GFLOPS, so I'm about 10% efficient with FPU
// This shows that I'm memory limited
//
// When using new code with macrofied arrays and TIMEORDER==2 and RK2 but still with PARALINE, I get 0.7GFLOPS since more focused on inversion than memory accesses.
//


////////////////////////////////////
//
// MPI Notes
//
////////////////////////////////////

// 1) Generally probably need to avoid fork() or system() calls
// 2) MPI-2 not fully supported on all systems, so ROMIO not fully supported

////////////////////////////////////
//
// OpenMP Notes
//
////////////////////////////////////

// 1) Ensure all called functions do not use global variables or else shared and data collision is possible.
//    This includes avoiding statically declared variables outside functions within a file.
// 2) omp construct can only be used on each part of loop, not on a multi-loop macro.
// 3) Use Intel(R) Thread Checker ?.? for Linux to check for OpenMP correctness.
//    NOT perfect: Missed "static FTYPE hslope" in coord.c when outside function
//
// a) source /opt/intel/itt/tcheck/bin/32e/tcheckvars.sh
// b) Compile code [for best results compile with -O -g
//    and set to run a finite number of steps]
// c) tcheck_cl ./grmhd
// Then look are report for data race conditions
//
// 4) On SMP systems, use less then total CPUs since basic system processes need some CPU time.  For example, with 1024 CPUs maybe use 1022.
//    On MPP or clusters, probably ok to use all systems since each system has same basic processes.
//
// 5) omp parallel construct should be far outside multi-dimen loops, with (e.g.) omp for construct on inner-most loop.
//
// 6) use the guided or dynamic schedule when each thread might take substantially different amounts of time, like when doing inversion.  Some threads might hit failure regions and take a long time.  Drawback is these have higher overheads (about 3X for chunk=10) compared to a static schedule.
//
// 7) Unlike in a single thread, across threads avoid cache coherence if updating nearby memory.  This is a problem known as false sharing.  If one thread updates part of memory in cache, that invalidates that cache line forcing other threads to reget the cache line from memory.
//
// 8) Note that if use pointer referencing outside parallel region and make that private (e.g. for *p2interp_l or *p2interp_r) inside the parallel region, then any prior pointer assignment is not preserved.  This is necessary since otherwise each thread would point to the same address.  So one must reassign that pointer relationship inside the parallel region.
//
// 9) For a 3D LOOP, if one wants to use "nowait", one has to use the "static" schedule so that upon leaving one loop and entering another, the thread is given the same exact seqeuence of iteration numbers.  Also assumes the loops (is,ie,js,je,ks,ke, etc.) from one construct to another are identical!  In either case, otherwise, the data set on a prior loop won't be available for the new loop for any given thread.
// So can only use nowait if both static and loop structures are identical.  This is rare.  And indeed, if static used across multiple loop blocks, then unlikely to benefit from nowait since static forces uniformity of workload.  So unlikely to be big difference in timing to reach end of prior loops.
//
// One can where *can* use nowait is when fully independent accesses.  This occurs, for example, when only 1 3D loop, but exterior loop is a DIMENLOOP/DIRLOOP.  In this case, each entrance to the next loop writes to totally different memory regions that depend upon dir itself (e.g. fluxvec(dir)).
//
// Again, for "nowait," must ensure that each loop does not use data written on previous loop.  If everything is independent, then "nowait" is ok.

// 10) Note that "static" variables inside functions are global in scope, so unless constant forever are not thread safe.
//     Typical "static int firsttime=1;  if(firsttime){ dosomething;firsttime=0;} cannot be used!
//     See metric_tools.c:matrix_inverse() for example of how deal with this.
// Can't have static in reconstructeno.c for memories (is then performance hit?)
// tau_neededbyharm.c from f2c: MUST remove static in front of variables!
//
//    Find these by doing: grep "  static" *.c *.h | less
//
//     Apart from static variables inside functions, static variables global to a single file can be dangerous if ever written to with multiple threads.  For example, coord.c has many static variables, but these should never be set when many threads are calling the coord.c's functions.
//
//    Find these by doing: grep "^static" *.c *.h | grep -v "(" | less
//
// Ensure that all uses of static are thread safe.
//
// 11) Must also ensure that "global variables" that aren't static (can be floating around inside files like ranc.c had "int called" as global, are either defined part of the global private or avoid global variable completely.
// Must check each file for such global varaibles.
// Use codenonfunc.sh and its extractnonfunc.c to help with this.
//
// 12) OpenMP Overhead:
//
//  Overhead on 64x32 calculation is about 22% with 78%=zcpsopenmp1core/zcps1core=38K/49K (TIMEORDER==2 FLUXB==FLUXCTTOTH).
// Saw no difference between static and guided schedules.
// No changes for static chunk size from 10 to 20 to 400.
// Currently operate at 36% efficiency with 4 cores, while 22% per core overhead would lead to a total of 78% efficiency.
//  So memory contention probably accounts for the other factor of 2-3X drop in efficiency.
//
//  Related to this, if number of things in loop is small, then 1 CPU can be faster than using OpenMP because of the overhead of a parallel region being introduced.  One can use the "if()" after "schedule" to avoid using the pragma if a certain condition is met.
//
// E.g. #pragma omp parallel for private(i__3,i__,j) schedule(guided) if(blocksize>100)
//
// As above, I noticed 1 core going with OpenMP is 22% slower than no-OpenMP going.  So one could eliminate the for loop pragma using:
// E.g. #pragma omp parallel for schedule(guided) if(numopenmpthreads>1)
// But we leave this out so that we can test OpenMP overhead since user is not expected to choose only 1 thread.
//
//
// 13) Volatile variable not allowed in parallel region since not thread safe.
//
//
// 14) private() and firstprivate() have overhead, so better to use scope to contain local variables.
// See: http://portal.acm.org/citation.cfm?doid=563647.563656
//      http://www.springerlink.com/content/680171t137073007/
//
// 15) copyin() can incur a LARGE fixed/core overhead if copyin() a large amount of data relative to work done.  Hit is 25%/core right now and that's pretty bad.

////////////////////////////////////
//
// Some related web resources for above comments
//
////////////////////////////////////


// http://www.kernel.org/pub/linux/kernel/v2.6/
// http://perfmon.sourceforge.net/
// http://sourceforge.net/projects/perfmon/
// http://perfmon2.sourceforge.net/
// http://perfmon2.sourceforge.net/
// http://perfsuite.ncsa.uiuc.edu/
// https://andrzejn.web.cern.ch/andrzejn/
// http://www.cyberciti.biz/tips/compiling-linux-kernel-26.html
// http://www.linuxhq.com/patch-howto.html
// https://andrzejn.web.cern.ch/andrzejn/#software
// http://www.eventhelix.com/RealtimeMantra/Basics/OptimizingCAndCPPCode.htm
// http://www.openwatcom.org/index.php/Performance_tuning
// http://www.netlib.org/eispack/
// http://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=2&t=178&start=0&st=0&sk=t&sd=a
// http://docs.hp.com/en/B3906-90006/index.html
// http://www.kde.org/
// http://sourceforge.net/projects/perfctr/
// http://perfmon2.sourceforge.net/pfmon_usersguide.html
// http://www.linuxjournal.com/article/7468
// http://brittany.nmsu.edu/~nmcac-intel/vtune/doc/users_guide/mergedProjects/analyzer_ec/mergedProjects/reference_olh/pentium4_hh/advice4_hh/dtlb_page_walk_misses.htm
// http://brittany.nmsu.edu/~nmcac-intel/vtune/doc/users_guide/mergedProjects/analyzer_ec/mergedProjects/reference_olh/pentium4_hh/advice4_hh/avoiding_dtlb_page_walk_misses.htm
// http://www.delorie.com/gnu/docs/binutils/gprof_25.html
// http://www.intel.com/software/products/mkl/data/vml/functions/pow.html
// http://www.eventhelix.com/RealtimeMantra/Basics/OptimizingCAndCPPCode.htm
// http://www.devx.com/SpecialReports/Door/40893
// http://www.ncsa.uiuc.edu/UserInfo/Resources/Software/Tools/PAPI/
// http://en.wikipedia.org/wiki/FLOPS
// http://www.maxxpi.net/pages/result-browser/top10---flops.php
// http://gcc.gnu.org/onlinedocs/cpp/Macros.html
// http://jamesthornton.com/emacs/node/emacs_97.html
// http://xahlee.org/emacs/find_replace_inter.html
// http://www.cs.utah.edu/dept/old/texinfo/emacs18/emacs_20.html
// http://www.emacswiki.org/emacs/SearchBuffers
// http://www.gnu.org/software/emacs/emacs-lisp-intro/html_node/emacs.html#Tags
// http://analyser.oli.tudelft.nl/regex/
// http://funarg.nfshost.com/r2/notes/sed-return-comma.html
// http://soft.zoneo.net/Linux/remove_empty_lines.php
// http://software.intel.com/en-us/articles/non-commercial-software-download/
// http://publib.boulder.ibm.com/infocenter/lnxpcomp/v7v91/index.jsp?topic=/com.ibm.vacpp7l.doc/language/ref/clrc09numnum.htm
// http://www.developers.net/intelisnshowcase/view/2180
// http://software.intel.com/en-us/intel-vtune/



// SVN comments:
//
// Initial HARM import:
// svn import harmjon http://harm.unfuddle.com/svn/harm_harm/ -m "Initial import"
// Additional directory import:
// svn import r8toras http://harm.unfuddle.com/svn/harm_harm/r8toras/ -m "r8toras import"
// Other imports:
// svn import sm2_4_1_jon http://harm.unfuddle.com/svn/harm_sm/ -m "Initial import"
// svn import vis5dplus_jon http://harm.unfuddle.com/svn/harm_vis5dplus/ -m "Initial import"
// svn import otherjoncodes http://harm.unfuddle.com/svn/harm_utilities/ -m "Initial import of other utilities"
//
//
// Then checkout to start svn tracking:
//
// svn checkout http://harm.unfuddle.com/svn/harm_harm/
//
//
// Added "r8toras" directory with files already inside:
//
// svn import r8toras http://harm.unfuddle.com/svn/harm_harm/r8toras/ -m "r8toras import"
// rm -rf r8toras
// svn update
//
//
// track user contributions:
// http://www.statsvn.org/
//
// diff without white space (kinda words?) 
// svn diff --diff-cmd diff -x -uw /path/to/file

// To revert a single file:
// svn merge -r <last revision>:<wanted revision> <filename>
// check: svn diff -r <wanted revision> <filename>



// git comments:
//

// Got git with:
// git@harm.unfuddle.com:harm/harmgit.git




// Setting up ssh tunnel (e.g. to get out of orange.slac.stanford.edu to portal through ki-rh42):
// On orange do:
// 1) ssh -N -f -L 2134:portal.astro.washington.edu:22 jon@ki-rh42.slac.stanford.edu  
// then still on orange do to ssh to jdexter@portal.astro.washington.edu through ki-rh42:
// 2) ssh -p 2134 jdexter@localhost
// OR to scp
// 3) scp -P 2134 dump???? jdexter@localhost:/astro/net/scratch1/jdexter/mb09/quadrupole_rout40/
// In either case you'll need to enter the jdexter password on portal.astro.washington.edu.



// Couldn't get this to work:
//
// Setting up ssh tunnel (e.g. to get data from queenbee to ranch through lonestar)
// On queenbee do (lonestar password):
// 1) ssh -N -f -L 2134:ranch.tacc.utexas.edu:22 tg802609@lonestar.tacc.utexas.edu
// then still on queenbee do to ssh to tg802609@ranch.tacc.utexas.edu through lonestar:
// 2) ssh -p 2134 tg802609@localhost
// OR to scp
// 3) scp -P 2134 dump???? tg802609@localhost:
// OR to bbcp
// 4) bbcp -k -a -r -P 5 -S 'ssh -p 2134 -x -a %I -l %U %H bbcp' dump???? tg802609@localhost:

// OR to bbcp w/ tarification
// 5) ??????????tar cvf - thickdisk3 | bbcp -P 5 -k -a -T 'ssh -x -a %I -l %U %H bbcp' tg802609@localhost:
// In latter cases you'll need to enter the tg802609 password on ranch


// This works:
//
// or directly:
// bbcp -r -k -a -P 5 thickdisk3_qb.tgz tg802609@ranch.tacc.utexas.edu:


// Scratch:
// tar cvf - thickdisk1 | ssh ${ARCHIVER} "cat > ${ARCHIVE}/thickdisk1.tar"
// ~/bin/bbcp -r -P 5 -a -k -T 'ssh -x -a %I -l %U %H bbcp' jon@ki-rh42.slac.stanford.edu:/media/disk/jon/thickdisk1/dumps/fieldline* jon@ki-rh42.slac.stanford.edu:/media/disk/jon/thickdisk1/dumps/dump0000* .

