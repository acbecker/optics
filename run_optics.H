# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <assert.h>
# include <iostream>
# include "misc.h"

/* HTM sizing */
typedef double		        float64;
#ifdef _WIN32
typedef __int64			int64;
typedef unsigned __int64	uint64;
typedef __int32                 int32;
typedef unsigned __int32        uint32;
#else
typedef long long		int64;
typedef unsigned long long	uint64;
typedef long                    int32;
typedef unsigned long           uint32;
#define IDHIGHBIT  0x8000000000000000LL
#define IDHIGHBIT2 0x4000000000000000LL
#endif

extern globalData GD;

/* -------------------------------------------------------------------------------- */
/*                                       lean.cc                                    */
/* -------------------------------------------------------------------------------- */
/*   (c) University of Munich, Database Group                                       */
/*   written by C. B"ohm, boehm@dbs.informatik.uni-muenchen.de                      */
/* -------------------------------------------------------------------------------- */
/* This program implements a very simplified version of the x-tree version published*/
/* in "Berchtold S., Keim D. A., Kriegel H.-P.: The X-Tree: An Index Structure for  */
/* High-Dimensional Data, Proc. 22th Int. Conf. on Very Large Data Bases, Bombay,   */
/* India, 1996, pp. 28-39."                                                         */
/* -------------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------------- */
/*                                       doublepq.h                                 */
/* -------------------------------------------------------------------------------- */
/*   (c) University of Munich, Database Group                                       */
/*   written by C. B"ohm, boehm@dbs.informatik.uni-muenchen.de                      */
/* -------------------------------------------------------------------------------- */
/* Helper functions for lean.cc                                                     */
/* -------------------------------------------------------------------------------- */

/*               doublepq.h                   */


template <class DTYPE> class FixDoublePQ {
    
private:
    int blocksize ;
    int num_elements ;
    
public:
    DTYPE * array ;
    int max_elements () {
        return (blocksize - sizeof (int) ) / sizeof (DTYPE) ;
    }
    
    FixDoublePQ () {
        blocksize = 65536 ;
        num_elements = 0 ;
        array = (DTYPE *) malloc (blocksize - sizeof (int) ) ;
    }
    
    FixDoublePQ (int _max_elements) {
        blocksize = sizeof (DTYPE) * _max_elements + sizeof (int) ;
        num_elements = 0 ;
        array = (DTYPE *) malloc (blocksize - sizeof (int) ) ;
    }
    
    ~FixDoublePQ () {
        free (array) ;
    }
    
    void write (FILE * fp, int blockno) {
# ifdef IO_MESSAGES
        fprintf (stderr, "w %2d %2d\n", blockno, (int) this % 100) ;
# endif
        fseek  (fp, blockno * blocksize, SEEK_SET) ;
        fwrite (& num_elements, sizeof (int), 1, fp) ;
        fwrite (array, blocksize - sizeof (int), 1, fp) ;
    }
    
    void read  (FILE * fp, int blockno) {
# ifdef IO_MESSAGES
        fprintf (stderr, "r %2d %2d\n", blockno, (int) this % 100) ;
# endif
        fseek  (fp, blockno * blocksize, SEEK_SET) ;
        fread  (& num_elements, sizeof (int), 1, fp) ;
        fread  (array, sizeof (DTYPE), max_elements (), fp) ;
    }
    
    void clear () {
# ifdef IO_MESSAGES
        fprintf (stderr, "c    %2d\n", (int) this % 100) ;
# endif
        num_elements = 0 ;
    }
    
    bool operator = (FixDoublePQ <DTYPE> & d) {
        blocksize = d.blocksize ;
        num_elements = d.num_elements ;
        //memcpy (array, d.array, num_elements * sizeof (DTYPE) ) ;
        /* ACB */
        return memcpy (array, d.array, num_elements * sizeof (DTYPE) ) ;
    }
    
    int get_num () {
        return num_elements ;
    }
    
    int empty () {
        return num_elements == 0 ;
    }
    
    int full () {
        return num_elements >= max_elements () ;
    }
    
    DTYPE * get_first () {
        if (num_elements == 0) {
            fprintf (stderr, "FixDoublePQ::get_first(): No element availlable\n") ;
            exit(-1);
        }
        return array ;
    }
    
    DTYPE * get_last () {
        if (num_elements == 0) {
            fprintf (stderr, "FixDoublePQ::get_last(): No element availlable\n") ;
            exit(-1);
        }
        if (num_elements == 1)
            return array ;
        return array + 1 ;
    }
    
    void dequeue_first (DTYPE * new_el) {
        DTYPE fill_el = new_el ? * new_el : array [-- num_elements] ;
        int i = 0 ;
        while (2*i+2 < num_elements &&
               (array [2*i+2] < fill_el || 
                (2*i+4 <= num_elements && array [2*i+4] < fill_el) ) ) {
            if (array [i+1] < fill_el) {
                array [i]   = array [i+1] ;
                array [i+1] = fill_el ;
                fill_el     = array [i] ;
            }
            array [i] = array [
                               i = num_elements < 2 * i + 4 || array [2*i+2] < array [2*i+4] 
                               ? 2*i+2 : 2*i+4] ;
        }
        if (i+1 <= num_elements && array [i+1] < fill_el) {
            array [i]   = array [i+1] ;
            array [i+1] = fill_el ;
        } else
            array [i] = fill_el ;
    }
    
    void dequeue_first (void) {
        dequeue_first (NULL) ;
    }
    
    void dequeue_last (DTYPE * new_el) {
        DTYPE * fill_el = new DTYPE (new_el ? * new_el : array [-- num_elements]) ;
        //	* fill_el = (new_el ? * new_el : array [-- num_elements]) ;
        int i = 0 ;
        while (2*i+3 < num_elements &&
               (* fill_el < array [2*i+3] || 
                (2*i+5 <= num_elements && * fill_el < array [2*i+5] ) ) ) {
            if (* fill_el < array [i]) {
                array [i+1] = array [i] ;
                array [i]   = * fill_el ;
                * fill_el     = array [i+1] ;
            }
            array [i+1] = array [1 + (
                                      i = num_elements < 2 * i + 5 || array [2*i+5] < array [2*i+3]
                                      ? 2*i+2 : 2*i+4) ] ;
        }
        if (2*i+2 < num_elements && * fill_el < array [2*i+2]) {
            array [i+1] = array [2*i+1] ;
            array [2*i+1] = * fill_el ;
        } else if (i+1 <= num_elements && * fill_el < array [i]) {
            array [i+1] = array [i] ;
            array [i] = * fill_el ;
        }
        else 
            array [i+1] = * fill_el ;
        delete fill_el ;
    }
    
    void dequeue_last (void) {
        dequeue_last (NULL) ;
    }
    
    void inqueue (DTYPE * new_el) {
        int i = num_elements ;
        if (num_elements == max_elements () ) {
            fprintf (stderr, "FixDoublePQ::inqueue(): Priority queue is full.\n") ;
            exit(-1);
        }
        num_elements ++ ;
        if (i % 2 && * new_el < array [i-1])
            array [i] = array [--i] ;
        while (i>1) {
            if (* new_el < array [(int)((i-2)/4)*2])
                array [i] = array [i=(int)((i-2)/4)*2] ;
            else if (array [(int)((i-2)/4)*2+1] < * new_el)
                array [i] = array [i=(int)((i-2)/4)*2+1] ;
            else
                break ;
        }
        array [i] = * new_el ;
    }
} ;

inline float64 fpow2 (float64 x) {
    register float64 h = x ;
    return h * h ;
}

class Data {
public:
    int dimension ;
    float64 * point ;
    /* ACB */
    uint64 id;
    
    Data (int d) {
        dimension = d ;
        point = new float64 [d] ;
        id = 0;
    }
    
    Data (int d, uint64 _id) {
        dimension = d ;
        point = new float64 [d] ;
        id = _id;
    }
    
    ~Data (void) {
        delete [] point ;
    }
    
    uint64 getId() const {
        return id;
    }
    void setId(uint64 _id) {
        id = _id;
    }
    
    /*
      double dist (double * q) {
      int i ;
      double h = 0 ;
      for (i=0 ; i<dimension ; i++)
      h += fpow2 (point [i] - q [i]) ;
      return sqrt (h) ;
      }
    */
} ;

template <class DTYPE> class DDist {
public:
    DTYPE * data ;
    float64 dist ;
    bool propagate ;
    
    //    DDist (void) {
    //	printf ("%%%d\n", (int) this) ;
    //    }
    
    DDist (DDist<DTYPE> * d) {
        //	printf ("^%d\n", (int) this) ;
        data = d->data ;
        dist = d->dist ;
        propagate = d->propagate ;
    }
    
    DDist (DTYPE * _data, float64 * q) {
        int i ;
        data = _data ;
        dist = data -> dist (q) ;
        propagate = true ;
        //	printf ("*%d\n", (int) this) ;
    }
    
    DDist (float64 _dist) {
        data = NULL ;
        dist = _dist ;
        propagate = false ;
        //	printf ("#%d\n", (int) this) ;
    }
    
    ~DDist () {
        //	printf ("/%d\n", (int) this) ;
    }
    
    bool operator < (DDist <DTYPE> & d) {
        return dist < d.dist ;
    }
} ;

template <class DTYPE> class DirEntry {
public:
    int dimension ;
    int num_entries ;
    float64 * lb ;
    float64 * ub ;
    DTYPE ** data ;
    
    DirEntry (DTYPE ** _data, int _num_entries) {
        int i, j ;
        
        data = _data ;
        dimension = (**data).dimension ;
        num_entries = _num_entries ;
        lb = new float64 [dimension] ;
        ub = new float64 [dimension] ;
        for (j=0 ; j<dimension ; j++)
            lb [j] = ub [j] = data [0] -> point [j] ;
        for (i=1 ; i<num_entries ; i++)
            for (j=0 ; j<dimension ; j++) {
                if (lb [j] > data [i] -> point [j])
                    lb [j] = data [i] -> point [j] ;
                if (ub [j] < data [i] -> point [j])
                    ub [j] = data [i] -> point [j] ;
            }
    }
    
    ~DirEntry (void) {
        delete lb ;
        delete ub ;
    }
    
    /*
      double dist (double * q) {
      int i ;
      double h = 0.0 ;
      for (i=0 ; i<dimension ; i++)
      if (q[i] < lb[i])
      h += fpow2 (lb[i] - q[i]) ;
      else if (q[i] > ub[i])
      h += fpow2 (q[i] - ub[i]) ;
      return sqrt (h) ;
      }
    */
    
    int maxextdim () {
        int i ;
        int m=0 ;
        for (i=1 ; i<dimension ; i++)
            if (ub [i] - lb [i] > ub [m] - lb [m])
                m = i ;
        return m ;
    }
    
    void out (FILE * outfile) {
        int i ;
        fwrite (& num_entries, sizeof (num_entries), 1, outfile) ;
        fwrite (lb, sizeof (float64), dimension, outfile) ;
        fwrite (ub, sizeof (float64), dimension, outfile) ;
        for (i=0 ; i<num_entries ; i++)
            fwrite (data [i] -> point, sizeof (float64), dimension, outfile) ;
    }
} ;

extern int sortdim;

inline int cmpdouble (float64 * f1, float64 * f2) {
    return *f1 > *f2 ? 1 : *f1 < *f2 ? -1 : 0 ;
}
inline int cmpdata (Data ** d1, Data ** d2) {
    float64 f1 = (**d1).point [sortdim] ;
    float64 f2 = (**d2).point [sortdim] ;
    return f1 > f2 ? 1 : f1 < f2 ? -1 : 0 ;
}

template <class DTYPE> class Leanindex {
public:
    int dimension ;
    int num_all ;
    int num_pages ;
    int server_id ;
    
    DTYPE ** dataset ;
    int dataset_m ;
    
    DirEntry <DTYPE> ** directory ;
    
    /* THIS SHIT WAS KILLING ME! */
    void wait_for_datapage (void) {}
    typedef int cmpfunc (const void *, const void *);
    /* THIS SHIT WAS KILLING ME! */
    
    Leanindex (int d, int _server_id) {
        dimension = d ;
        server_id = _server_id ;
        num_all = 0 ;
        num_pages = 0 ;
        directory = 0 ;
        dataset_m = 10000 * sizeof (DTYPE *) ;
        dataset = (DTYPE **) malloc (dataset_m) ;
        if (! dataset) {
            printf ("Error in malloc. Exiting\n\n") ;
            exit (-1) ;
        }
    }
    
    ~Leanindex (void) {
        if (dataset) free (dataset) ;
        if (directory) free (directory) ;
    }
    
    void out (FILE * outfile) {
        int i ;
        int maxc = 0 ;
        if (! directory) {
            fprintf (stderr, "Fehler, kein Directory\n\n") ;
            exit (-1) ;
        }
        for (i=0 ; i<num_pages ; i++)
            if (maxc < directory[i]->num_entries)
                maxc = directory[i]->num_entries ;
        fwrite (& dimension, sizeof (int), 1, outfile) ;
        fwrite (& maxc, sizeof (int), 1, outfile) ;
        fwrite (& num_all, sizeof (int), 1, outfile) ;
        fwrite (& num_pages, sizeof (int), 1, outfile) ;
        for (i=0 ; i<num_pages ; i++)
            directory[i]->out (outfile) ;
    }
    
    void insert (DTYPE * data) {
        if (directory) {
            printf ("Warning: No dynamic insert after construct.\n"
                    "Discarding Directory\n\n") ;
            free (directory) ;
            directory = 0 ;
        }
        if (! dataset) {
            printf ("Error: No dataset allocated. Exiting\n\n") ;
            exit (-1) ;
        }
        if ((num_all+1) * sizeof (DTYPE) > dataset_m) {
            dataset_m = dataset_m * 3 / 2 ;
            dataset = (DTYPE **) realloc (dataset, dataset_m) ;
            if (! dataset) {
                printf ("Error in malloc. Exiting\n\n") ;
                exit (-1) ;
            }
        }
        dataset [num_all] = data ;
        num_all ++ ;
    }
    
    typedef struct {
        float64 extension ;
        int dimno ;
        float64 lb ;
        float64 ub ;
    } diminfo ;
    
    void construct_rec (int capacity, int start, int end) {
        int i, j, k ;
        DirEntry <DTYPE> * de = new DirEntry <DTYPE> (dataset + start, end - start) ;
        DTYPE * h ;
        if (end - start <= capacity) {
            directory [num_pages] = de ;
            num_pages ++ ;
            diminfo * di = (diminfo *) malloc (sizeof (diminfo) * dimension) ;
            for (i=0 ; i<dimension ; i++) {
                di [i] . dimno     = i ;
                di [i] . lb        = de -> lb [i] ;
                di [i] . ub        = de -> ub [i] ;
                di [i] . extension = de -> ub [i] - de -> lb [i] ;
            }
            qsort (di, dimension, sizeof (diminfo), (cmpfunc *) cmpdouble) ;
            k = 0 ;
            for (j=0 ; j<dimension ; j++) {
                for (i=k ; i<de -> num_entries ; i++)
                    if (di [j] . lb == dataset[start+i]->point[di[j].dimno]) {
                        h = dataset[start+i] ;
                        dataset[start+i] = dataset[start+k] ;
                        dataset[start+k] = h ;
                        k ++ ;
                        break ;
                    }
                for (i=k ; i<de -> num_entries ; i++)
                    if (di [j] . ub == dataset[start+i]->point[di[j].dimno]) {
                        h = dataset[start+i] ;
                        dataset[start+i] = dataset[start+k] ;
                        dataset[start+k] = h ;
                        k ++ ;
                        break ;
                    }
            }
        } else {
            sortdim = de -> maxextdim () ;
            delete de ;
            qsort (dataset + start, end - start, sizeof (DTYPE *), (cmpfunc *) cmpdata) ;
            i = (end - start) / capacity + ((end - start) % capacity ? 1 : 0) ;
            j = (int) ((float64) (end - start) * (float64) (int) (i / 2) / (float64) i + 0.5);
            construct_rec (capacity, start, start + j) ;
            construct_rec (capacity, start + j, end) ;
        }
    }
    
    void construct (int capacity) {
        int i = num_all / capacity + (num_all % capacity ? 1 : 0) ;
        directory = new DirEntry <DTYPE> * [i] ;
        construct_rec (capacity, 0, num_all) ;
    }
    
    // some communication routines:
    
    void get_distances (FixDoublePQ <DDist <DTYPE> > * candidates) {
        int ne = 0 ; // Number of Entries received
        int i ;
        DDist <DTYPE> * dd ;
        /***********
	for (i=0 ; i<ne ; i++)
	    if (source [i] != server_id) {
	        dd = new DDist <DTYPE> (dist [i]) ;
		if (! candidates.full () )
		    candidates.inqueue (dd) ;
		else if (dd -> dist < candidates.get_last () -> dist)
		    candidates.dequeue_last (dd) ;
		delete (dd) ;
            }
        ************/
    }
    
    void put_distances (FixDoublePQ <DDist <DTYPE> > * candidates) {
        int ne = 0 ; // Number of Entries
        int i ;
        for (i = 0 ; i < candidates -> get_num () ; i ++)
            if (candidates -> array [i] . propagate) {
                // write server_id
                // write candidates -> array [i] . dist
                ne ++ ;
                candidates -> array [i] .propagate = false ;
            }
    }
    
    void put_resume () {
        // write server_id
    }
    
    void knn_search_hs (float64 * q, int k) {
        FixDoublePQ <DDist <DirEntry <DTYPE> > > apl (num_pages) ;
        FixDoublePQ <DDist <DTYPE> > candidates (k) ;
        DDist <DirEntry <DTYPE> > * dde ;
        DirEntry <DTYPE> * de ;
        DDist <DTYPE> * ddata ;
        DTYPE * data ;
        int i, page_access=0 ;
        
        // Fill the Active Page List (apl) with MINDIST values for all Data Page Regions
        for (i=0 ; i<num_pages ; i++) {
            dde = new DDist <DirEntry <DTYPE> > (directory [i], q) ;
            apl.inqueue (dde) ;
            delete dde ;
        }
        
        // Main Loop of the HS algorithm
        while (! apl.empty () 
               && (! candidates.full () 
                   || apl.get_first () -> dist < candidates.get_last() -> dist) ) {
            // Process Data Page apl.get_first
            
            wait_for_datapage () ;
            page_access++;
            de = apl.get_first() -> data ;
            for (i = 0 ; i < de->num_entries ; i ++) {
                ddata = new DDist <DTYPE> (de->data[i], q) ;
                //		printf ("+%d \n", (int) ddata) ;
                if (! candidates.full () )
                    candidates.inqueue (ddata) ;
                else if (ddata -> dist < candidates.get_last () -> dist) {
                    candidates.dequeue_last () ;
                    candidates.inqueue (ddata) ;
                }
                //		printf ("-%d \n", (int) ddata) ;
                delete ddata ;
            }
            apl.dequeue_first () ;
            // Exchange Distance Information with Coordinator
            get_distances (& candidates) ;
            put_distances (& candidates) ;
        }
        // Report the Results
        put_resume () ;
        /*
          for (i=0 ; i<candidates.get_num () ; i++)
          printf ("%4.4f, ", candidates.array [i]. dist) ;
          printf ("\n") ;
        */
        while (! candidates.empty () ) {
            printf ("%4.4f, ", candidates.get_first () -> dist) ;
            candidates.dequeue_first () ;
        }
        printf ("\n") ;
        printf("page_access = %d\n", page_access);
    }
} ;


