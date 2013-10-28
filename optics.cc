/* -------------------------------------------------------------------------------- */
/*                                       OPTICS                                     */
/* -------------------------------------------------------------------------------- */
/*   (c) University of Munich, Database Group                                       */
/*   written by M. Breunig, breunig@dbs.informatik.uni-muenchen.de                  */
/* -------------------------------------------------------------------------------- */
/* This program implements a very simple version of the OPTICS algorithm published  */
/* in "Ankerst M., Breunig M. M., Kriegel H.-P., Sander J.: OPTICS: Ordering Points */
/* To Identify the Clustering Structure, Proc. ACM SIGMOD Int. Conf. on Management  */
/* of Data, Philadelphia, PA, 1999."                                                */
/* -------------------------------------------------------------------------------- */

#include <iostream>
#include <iomanip.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include "optics.H"

// --------------------------------------------------------------------------------
// execute the optics algorithm
optics::optics(dataBase *_db) {
    int donePoints = 0;
    int index;
    float64 r, coreLevel;
    ptCluster *pc = new ptCluster;
    
    nclusters = 0;
    //clusters = (ptCluster **)malloc(sizeof(ptCluster *));
    //ptCluster clusters[100];
    
    db = _db;
    sl = new seedList(db->getNumPoints()); assert(sl);
    ndist = new float64[GD.minPts]; assert(ndist);
    while(donePoints<db->getNumPoints()) {
        if(sl->getNumEntries()==0) {
            // pick a random point and add to seedlist with r=infinity
            sl->addRandomPoint(db);
        }
        sl->getMin(&index, &r); // get and delete point with smallest r from seedlist
        coreLevel = updateAllReach(index); // update seedlist with reachabilities from point index
        //output(index, coreLevel, r); // output this point with it's current reachability
        ++donePoints;
        
        /* CAW */
        point p = db->getPoint(index);
        if (r == OPTICS_INF) {
            if (nclusters == 0) {
                /* first cluster, do nothing yet */
            }
            else {
                /* assign old one */
                clusters.push_back(pc);
            }	 
            nclusters += 1;
            pc = new ptCluster;
            pc->add(p, index);
        }
        else {
            pc->add(p, index);
        }
    }   
    /* for last cluster... */
    clusters.push_back(pc);
}

// --------------------------------------------------------------------------------
// free memory
optics::~optics() {
    delete sl;
    delete[] ndist;
}

// --------------------------------------------------------------------------------
// enter/update all the reachabilities of all points as seen from point index
// standard version of OPTICS (non-join based)
float64 optics::updateAllReach(int index) {
    ptListContent *pc;
    int i;
    
    for(i=0; i<GD.minPts; ++i) 
        ndist[i] = OPTICS_INF;
    
    // execute range query
    ptList *pl = db->rangeQuery(db->getPoint(index));
    
    // compute core level
    ptListIter ptli(pl);
    while((pc = ptli.next())) {
        if(ndist[GD.minPts-1]>pc->d) {
            // have to include this dist into ndist array
            i = GD.minPts-2;
            while((ndist[i]>pc->d)&&(i>=0)) {
                ndist[i+1] = ndist[i];
                --i;
            }
            ndist[i+1] = pc->d;
        }
    }
    float64 coreLevel = ndist[GD.minPts-1];
    
    // update all points in eps-neighborhood
    if(coreLevel<OPTICS_INF) {
        // point is core point
        ptListIter ptli(pl);
        while((pc = ptli.next())) {
            sl->update(pc->pts, max(pc->d, coreLevel));
        }
    }
    delete pl;
    return(coreLevel);
}

// --------------------------------------------------------------------------------
// output the point index with reachability r 
//void optics::output(int index, float64 c, float64 r) {
//    // cannot use stream functions to print floats, because
//    // even if we set ios::fixed, it sometimes prints e-notation.
//    const point &p = db->getPoint(index);
//    
//   for(int i=0; i<GD.dim; ++i) {
//        fprintf(stderr, "%10.8f ", p.getCoord(i));
//    }
//    fprintf(stderr, "%10.8f %10.8f \"%d\"", c, ((r==OPTICS_INF)?-1:r), index);
//    fprintf(stderr, "\n");
//}

// ================================================================================
// ========== class seedList                                             ==========
// ================================================================================

// init all data structures
seedList::seedList(int _size) {
    size       = _size;
    numEntries = 0;
    doneIndex  = 0;
    slReach    = new float64[size];   assert(slReach);
    slIndex    = new int[size];     assert(slIndex);
    slInvIndex = new int[size];     assert(slInvIndex);
    for(int i=0; i<size; ++i)
        slIndex[i] = NOTHANDLED;
}

seedList::~seedList() {
    delete[] slReach;
    delete[] slIndex;
    delete[] slInvIndex;
}

void seedList::addRandomPoint(dataBase *db) {
    // do NN-query around origin
    int index = db->specialNNQuery(slIndex);
    update(index, OPTICS_INF);
}

// adjust the heap structure, assuming only position "index" IN THE HEAP has been
// updated/added
void seedList::adjustHeap(int index) {
    int j;
    float64 hf;
    int hi;
    
    // do we need to bubble up or sink down?
    if(index>0) { // we are not the root
        j = (index-1)/2; // j is father
        if(slReach[index]<slReach[j]) {
            // bubble up
            while(index>0) {
                // exchange(index, j);
                hf = slReach[index];
                slReach[index] = slReach[j];
                slReach[j] = hf;
                
                slIndex[slInvIndex[index]] = j;
                slIndex[slInvIndex[j]] = index;
                
                hi = slInvIndex[index];
                slInvIndex[index] = slInvIndex[j];
                slInvIndex[j] = hi;
                // end exchange(index, j);
                index = j;
                if(index>0) {
                    j = (index-1)/2;
                    if(slReach[index]>=slReach[j]) // heap ok
                        index = 0;
                }
            }
            return; // STOP HERE
        }
    }
    // sink down if necessary
    j = 2*index+1; // left son
    while(j < numEntries) { // we are not a leaf
        if(j<numEntries-1)  // right leaf does exist
            if(slReach[j+1]<slReach[j]) // right leaf is smaller
                ++j;
        // slReach[j] is smaller son if slReach[index]!
        if(slReach[j]<slReach[index]) { // heap not ok
            // exchange(index, j);
            hf = slReach[index];
            slReach[index] = slReach[j];
            slReach[j] = hf;
            
            slIndex[slInvIndex[index]] = j;
            slIndex[slInvIndex[j]] = index;
            
            hi = slInvIndex[index];
            slInvIndex[index] = slInvIndex[j];
            slInvIndex[j] = hi;
            // end exchange(index, j);
            index = j;
            j = index*2+1;
        } else  // heap ok
            j = numEntries;
    }
}

// add or update point "index" with reachability value "r"
void seedList::update(int index, float64 r) {
    if(slIndex[index]!=HANDLED) {
        if(slIndex[index]==NOTHANDLED) {
            // point needs to be added to the seedlist
            slIndex[index] = numEntries;
            slReach[numEntries] = r;
            slInvIndex[numEntries] = index;
            numEntries++;
            adjustHeap(numEntries-1);
        } else {
            // point is already in the seedlist
            if(slReach[slIndex[index]] > r) {
                // need to update because new reachability is smaller
                slReach[slIndex[index]] = r;
                adjustHeap(slIndex[index]);
            }
        }
    }
}

// get and delete the point with the smallest r-value from the seedlist
void seedList::getMin(int *index, float64 *r) {
    *r = slReach[0];
    *index = slInvIndex[0];
    --numEntries;
    slIndex[slInvIndex[0]] = HANDLED;
    if(numEntries>0) {
        // do this only if the list is not empty
        slReach[0] = slReach[numEntries];
        slInvIndex[0] = slInvIndex[numEntries];
        slIndex[slInvIndex[0]] = 0;
        adjustHeap(0);
    }
}


// ================================================================================
// ========== class leanTreeDataBase                                             ==========
// ================================================================================

// --------------------------------------------------------------------------------
// prepare database to be read into memory
leanTreeDataBase::leanTreeDataBase(void) {
    page = NULL;
}

// --------------------------------------------------------------------------------
// if a database was read into memory, free the memory
leanTreeDataBase::~leanTreeDataBase() {
    if(page!=NULL) {
        delete[] page;
    }
}

// --------------------------------------------------------------------------------
// load the lean tree data in from the Leanindex structure, this marries leantree to optics!
int leanTreeDataBase::load(Leanindex <Data> leanindex) {
    int i,j, k, nRead=0, numEntries;
    point *p, *tp;
    
    printf("\n**** Optics Reading Ltree\n");
    GD.dim = leanindex.dimension;
    maxEntriesPerPage = 0;
    for (i=0; i < leanindex.num_pages; i++)
        if (maxEntriesPerPage < leanindex.directory[i]->num_entries)
            maxEntriesPerPage = leanindex.directory[i]->num_entries;
    
    numPoints = leanindex.num_all;
    numPages = leanindex.num_pages;
    
    p = new point;
    page = new dataPage[numPages];                             assert(page);
    // read the data pages
    for(i=0; i<numPages; ++i) {
        tp = new point;
        numEntries = leanindex.directory[i]->num_entries;
        //fprintf(stderr, "DEBUG : %d %d\n", i, numEntries);
        page[i].setNumEntries(numEntries);
        page[i].firstPointIndex = nRead;
        
        for (k=0; k<GD.dim; k++) {
            tp->setCoord(k, leanindex.directory[i]->lb[k]);
        }
        page[i].rect.setLB(*tp); 
        
        for (k=0; k<GD.dim; k++) {
            tp->setCoord(k, leanindex.directory[i]->ub[k]);
        }
        page[i].rect.setUB(*tp); 
        
        delete tp;
        
        for(j=0; j<numEntries; ++j) {
            tp = new point;
            for (k=0; k<GD.dim; k++) {
                tp->setCoord(k, leanindex.directory[i]->data[j]->point[k]);
            }
            tp->setId(leanindex.directory[i]->data[j]->getId());
            
            /*
              fprintf(stderr, "CAw %f %f %llu\n", 
              leanindex.directory[i]->data[j]->point[0],
              leanindex.directory[i]->data[j]->point[1],
              leanindex.directory[i]->data[j]->id);
            */
            
            page[i].entry[j].setPoint(*tp);
            page[i].rect.include(*tp);
            delete tp;
            ++nRead;
        }
    }
    printf("\r     numPagesRead      = %d\n",numPages);
    printf("     numDatapointsRead = %d\n",nRead);
    delete p;
    return(0);
}

// --------------------------------------------------------------------------------
// execute a range query on the database
// NOTE: the returned ptList has to be deleted by the caller!
ptList * leanTreeDataBase::rangeQuery(const point &p) {
    ptList *pl = new ptList;
    ptListContent pc;
    
    for(int i=0; i<numPages; ++i) { // for all pages
        if(page[i].rect.dist(p)<=GD.eps) { // page may contain points in range
            for(int j=0; j<page[i].numEntries; ++j) {
                if((pc.d=p.distSC(page[i].entry[j].getPoint()))<OPTICS_INF) {
                    pc.pts = page[i].firstPointIndex+j;
                    pl->add(&pc);
                }
            }
        }
    }
    return(pl);
}

// --------------------------------------------------------------------------------
// do a special NN-query: compute point next to origin that has not been
// handled
int leanTreeDataBase::specialNNQuery(int *handled) {
    float64 currdist=0, newdist;
    int index = -1, i, j;
    point origin;
    
    /* ACB */
    if ((GD.dtype == 1) && (GD.dim == 3)) {
        origin.setCoord(0, 1);
        for(i=1; i<GD.dim; ++i) origin.setCoord(i, 0);
    }
    else {
        for(i=0; i<GD.dim; ++i) origin.setCoord(i, 0);
    }
    
    
    for(i=0; i<numPages; ++i) { 
        // for every page
        if((page[i].rect.dist(origin) < currdist)||(index==-1)) {
            // consider page
            for(j=0; j<page[i].numEntries; ++j) {
                // for every point in the page
                if(handled[page[i].firstPointIndex+j]==NOTHANDLED) {
                    // point has not already been handled
                    newdist = origin.dist(page[i].entry[j].getPoint());
                    if((newdist<currdist)||(index==-1)) {
                        currdist = newdist;
                        index = page[i].firstPointIndex + j;
                    }
                }
            }
        }
    }
    return(index);
}


// ================================================================================
// ========== end of file                                                ==========
// ================================================================================
