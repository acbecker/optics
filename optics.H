#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <iostream>
#include <vector.h>
#include "run_optics.H"
#include "misc.h"

extern globalData GD;

#define assignMin(a,b) if(a>b) a=b;
#define assignMax(a,b) if(a<b) a=b;
#define UNDEF -1
#define TRUE (1==1)
#define FALSE (1==0)
#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))
#define OPTICS_INF (2*GD.eps)


// ======================================== 
// a data point
/* ACB - the way this is used, dist is never called when this is not point 0,0 */
/* this screws up the whole 3-d sperical coordinate system.  */
/* maybe need to make this 1, 0, 0 in that case... */
class point {
public:
    point(void) { 
        init(); 
    }
    point(const point &p) { 
        init(); assign(p); 
    }
    ~point() { 
        delete[] f; 
    }
    uint64 getId() const {
        return id;
    }
    void setId(uint64 _id) {
        id = _id;
    }
    float64 getCoord(int _dim) const { 
        //if ( ((id == 3) or (id == 5)) and _dim == 0) {
        //printE();
        //}

        assert((_dim>=0)&&(_dim<GD.dim)); 
        return(f[_dim]); 
    }
    void  setCoord(int _dim, float64 _f) { 
        assert((_dim>=0)&&(_dim<GD.dim)); 
        f[_dim]  = _f; 
    }
    void  assign(const point &p) { 
        for(int i=0; i<GD.dim; ++i) {
            f[i]  = p.getCoord(i); 
        }
        setId(p.getId());
    }
    
    void print(void) const { 
        std::cout << "[ "; 
        for(int i=0; i<GD.dim-1; ++i) {
            std::cout << f[i] << ", ";
        }
        std::cout << f[GD.dim-1] << "]\n" << std::flush;
    }

    void printE(void) const { 
        fprintf(stdout, "[ ");
        for(int i=0; i<GD.dim; ++i) {
            fprintf(stdout, "%.11e ", f[i]);
        }
        fprintf(stdout, "]\n");
    }

    void fileread(FILE *infile) { 
        fread(f, sizeof(float64), GD.dim, infile); 
        assert(!feof(infile)); 
    }
    /* ************************************************* */
    /* ACB - change from euclidean distance to spherical */
    /* ************************************************* */
    float64 dist(const point &p) const {
        if (GD.dtype == 0) {
            return distEuc(p);      
        }
        else if (GD.dtype == 1) {
            return distEuc(p);
            /* This was changed to use the normal Euclidiean metric */
            //return distCos(p);
        }
        else
            return GD.eps;
    }
    float64 distEuc(const point &p) const { 
        /* Euclidian Distance */
        float64 dist = 0., offset;
        for(int i=0; i<GD.dim; ++i) {
            offset = p.f[i]-f[i];
            dist  += (offset*offset);
        }
        return(sqrt(dist));
    }
    float64 distCos(const point &p) const{
        /* 3-D Euclidian Distance where distance metric is 1-cos(epsilon) */
        /* add you idiot, cos(epsilon) = a.dot.b, 
           meaning you multiply the terms together, not add! */
        float64 dist = 0., offset;
        for(int i=0; i<GD.dim; ++i) {
            offset = p.f[i]*f[i];
            dist  += offset;
            //fprintf(stdout, "#CAW4 %d %.11e %.11e %.11e %.11e %.11e\n", i, GD.eps, p.f[i], f[i], offset, dist);
        }
        //fprintf(stdout, "#CAW5 %e %e %d\n", dist, 1-dist, ((1-dist) < GD.eps) ? 1 : 0);
        return 1.-dist;
    }
    
    /* ************************************************* */
    /* ACB - change from euclidean distance to spherical */
    /* ************************************************* */
    float64 distSC(const point &p) const{
        if (GD.dtype == 0) {
            return distSC_Euc(p);
        }
        else if (GD.dtype == 1) {
            return distSC_Euc(p);
            /* This was changed to use the normal Euclidiean metric */
            //return distSC_Cos(p);
        }
        else return GD.eps;
    }
    float64 distSC_Euc(const point &p) const { 
        // distSC(p1,p2) if distance between p1 and p2 is smaller 
        // than or equal to GD.eps and OPTICS_INF otherwise
        float64 offset, eps;
        eps = GD.eps*GD.eps;
        for(int i=0; i<GD.dim; ++i) {
            offset = p.f[i]-f[i];
            eps   -= (offset*offset);
            if(eps<0) return(OPTICS_INF);
        }
        return(sqrt(GD.eps*GD.eps-eps));
    }
    float64 distSC_Cos(const point &p) const {
        /* take out what they did above for speed */
        float64 dist;
        dist = distCos(p);
        if (dist < GD.eps)
            return dist;
        return OPTICS_INF;
    }
    
private:
    inline void init(void) { 
        f = new float64[GD.dim]; 
    }
    float64  *f;     // values
    uint64 id;      // ID
    friend class mbr;
};

// ======================================== 
// a mbr (minimum bounding rectangle)
class mbr {
public:
    mbr(void) 
    { }
    mbr(const point &p) 
    { set(p); }
    ~mbr()
    { }
    void set(const point &p) 
    { setLB(p); setUB(p); }
    void setLB(const point &p) 
    { lb.assign(p); }
    void setUB(const point &p) 
    { ub.assign(p); }
    const point &getLB(void)
    { return(lb); }
    const point &getUB(void)
    { return(ub); }
    void include(const point &p) // resize the mbr so that it includes p
    { for(int i=0; i<GD.dim; ++i) {
            assignMin(lb.f[i], p.f[i]); 
            assignMax(ub.f[i], p.f[i]); 
        } 
    }
    /* ************************************************* */
    /* ACB - change from euclidean distance to spherical */
    /* ************************************************* */
    float64 dist(const point &p) {
        if (GD.dtype == 0) {
            return distEuc(p);
        }
        else if (GD.dtype == 1) {
            return distEuc(p);
            /* This was changed to use the normal Euclidiean metric */
            //return distCos(p);
        }
        else return GD.eps;
    }
    float64 distEuc(const point &p) { 
        float64 dist = 0, offset;
        for(int i=0; i<GD.dim; ++i) {
            if(ub.f[i] <= p.f[i])
                offset = p.f[i] - ub.f[i];
            else if(lb.f[i] >= p.f[i])
                offset = lb.f[i] - p.f[i];
            else
                offset = 0;
            dist += (offset*offset);
        }
        return(sqrt(dist));
    }
    /* ACB - again, its a dot product man! */
    float64 distCos(const point &p) {
        float64 dist = 0, offset;
        for(int i=0; i<GD.dim; ++i) {
            if(ub.f[i] <= p.f[i])
                offset = p.f[i] * ub.f[i];
            else if(lb.f[i] >= p.f[i])
                offset = lb.f[i] * p.f[i];
            else
                offset = 0;
            //offset = p.f[i]; /* ACB - why did I have this in there? */
            dist += (offset*offset);
        }
        return 1.-dist;
    }
    
private:
    point lb; // lower bound
    point ub; // upper bound
};

// ======================================== 
// a point list
#define PTLISTLEN 20
class ptListContent {
public:
    int pts; // distance of point pts 
    float64 d; // is d
};

class ptListEntry {
public:
    ptListEntry(void) { next = NULL; }
    ~ptListEntry()    { if(next) delete next; }
private:
    friend class ptList;
    friend class ptListIter;
    ptListContent pc[PTLISTLEN];
    ptListEntry *next;
};

class ptList {
public:
    ptList(void)   { nextfree = 0; entry = NULL; }
    ~ptList()      { if(entry) delete entry; }
    void add(ptListContent *pc) { 
        if(entry==NULL) { // list is completely empty
            entry = new ptListEntry;
        } else if(nextfree==PTLISTLEN) { // current entry is full
            ptListEntry *newentry = new ptListEntry;
            newentry->next = entry;
            entry = newentry;
            nextfree = 0;
        }
        entry->pc[nextfree] = *pc;
        ++nextfree;
    }
    
private:
    friend class ptListIter;
    int nextfree;
    ptListEntry *entry;
};

class ptListIter {
public:
    ptListIter(const ptList *pl) { 
        maxEntry = pl->nextfree; entry = pl->entry; i=0; 
    }
    ptListContent *next(void) { 
        if(entry==NULL) return(NULL);
        if(i<maxEntry) {
            ++i;
            return(&(entry->pc[i-1]));
        } else {
            entry = entry->next;
            i=0;
            maxEntry = PTLISTLEN;
            return(next());
        }
    }
private:
    int maxEntry, i;
    ptListEntry *entry;
};

class ptCluster {
public:
    vector<const point*> points;   
    vector<int> inds;
    float64 ra;
    float64 dec;
    int    npts;
    
    ptCluster(void) {
        ra = 0;
        dec = 0;
        npts = 0;
    }
    
    ~ptCluster() {
    }
    
    void add (const point &p, int index) {
        /* somehow this info gets lost by the time i get to center */
        //points.push_back(&p);
        
        /* inefficient but it works for god's sake */
        points.push_back(new point(p));
        inds.push_back(index);
        npts += 1;
        
        /* ACB : this works */
        //cout << "caw0 " << index << " " << p.getCoord(0) << " " << points.back()->getCoord(0) << endl;
        //cout << "caw1 " << index << " " << p.getCoord(1) << " " << points.back()->getCoord(1) << endl;
        //cout << "caw2 " << index << " " << p.getCoord(2) << " " << points.back()->getCoord(2) << endl;
    }
    
    void center() {
        int i;
        float64 dras=0, ddecs=0, tmpdec=0, tmpra=0;
        //const double rad2deg = 57.2957795130823208767981548141052;
        
        for (i=0; i<npts; i++) {
            //const point *p = points[i];
            
            /* all the pointers in points are empty by now...   sniff... */
            if (GD.dim == 2) {
                dras  += points.at(i)->getCoord(0);
                ddecs += points.at(i)->getCoord(1);
                //cout << "caw " << i << " " << points.at(i)->getCoord(0) << " " << points.at(i)->getCoord(1) << endl;
                //cout << "caw " << i << " " << points[i]->getCoord(0) << " " << points[i]->getCoord(1) << endl;
            }
            else if (GD.dim == 3) {
                // what the hell, this does not work!?
                // what the hell, this does not work!?
                // what the hell, this does not work!?

                tmpdec = asin( points.at(i)->getCoord(2));
                tmpra  = atan2(points.at(i)->getCoord(1), points.at(i)->getCoord(0));

                /* bad if cos = 0 */
                //dras  += asin(points.at(i)->getCoord(1) / cos(tmpdec));

                //dras  += atan2(points.at(i)->getCoord(1), points.at(i)->getCoord(0));
                //cout << " caw " << i << " " << points.at(i)->getCoord(0) ;
                //cout << " " << points.at(i)->getCoord(1);
                //cout << " " << points.at(i)->getCoord(2) << endl;
                //cout << " " << tmpdec * rad2deg;
                //cout << " " << atan2(points.at(i)->getCoord(1), points.at(i)->getCoord(0))*rad2deg;
                //cout << " " << acos(points.at(i)->getCoord(0) / cos(tmpdec))*rad2deg;
                //cout << " " << asin(points.at(i)->getCoord(1) / cos(tmpdec))*rad2deg << endl;

                ddecs += tmpdec;
                dras  += tmpra;
            }
        }
        
        ra  = dras / npts;
        dec = ddecs / npts;
    }
    
private:
    int NS;
    int nalloc;
};

// ======================================== 
// an entry in a datapage
class dataBaseEntry {
public:
    dataBaseEntry()              { }
    ~dataBaseEntry()             { }
    void setPoint(point &_p)     { p.assign(_p);  }
    const point &getPoint(void)  { return(p); }
private:
    point p;        // the data point
};

// ======================================== 
// a datapage
class dataPage {
public:
    dataPage(void)
    { numEntries = 0; entry = NULL; }
    ~dataPage()
    { if(entry!=NULL) delete[] entry; }
    void setNumEntries(int _numEntries)
    { numEntries = _numEntries;
        if(entry!=NULL) delete[] entry;
        entry = new dataBaseEntry[numEntries]; 
    }
    int             firstPointIndex;   // number of the first point in this page
    int             numEntries;        // how many entries does this data page have
    dataBaseEntry*  entry;             // array of size "numEntries"
    mbr             rect;              // minimal bounding rectangle of this page
private:
};

// ======================================== 
// the database - abstract base class
class dataBase {
public:
    dataBase(void) { }
    virtual ~dataBase() { }
    int                      getNumPoints(void) { return numPoints; }
    virtual const point &    getPoint(int index) = 0;
    virtual ptList *         rangeQuery(const point &p) = 0;
    virtual int              specialNNQuery(int *handled) = 0;
    //globalData GD;
protected:
    int numPoints;
};

// ======================================== 
// a very simple database, no index, just the points
class simpleDataBase : public dataBase {
public:
    simpleDataBase(void);
    ~simpleDataBase();
    int              load(float64 *x, float64 *y, int _numPoints);
    const point &    getPoint(int index) { assert(index<numPoints); return(entry[index].getPoint()); }
    ptList *         rangeQuery(const point &p);
    int              specialNNQuery(int *handled);
private:
    dataBaseEntry *entry;
};

// ======================================== 
// a database using a lean-tree
class leanTreeDataBase : public dataBase {
public:
    leanTreeDataBase(void);
    ~leanTreeDataBase();
    int load(Leanindex <Data> leanindex);
    const point &getPoint(int index) { 
        int p = pageContainingPoint(index);
        return(page[p].entry[index-page[p].firstPointIndex].getPoint());
    }
    ptList *         rangeQuery(const point &p);
    int              specialNNQuery(int *handled);
    int              getNumPages(void) { return numPages; }
    int              getPageSize(void) { return pageSize; }
    int              getMaxEntriesPerPage(void) { return maxEntriesPerPage; }
private:
    int pageContainingPoint(int index) {
        // compute pageNum containing point index
        assert(index<numPoints);
        int p=0, s=0;
        while(s+page[p].numEntries-1<index) { s += page[p].numEntries; ++p; }
        return(p);
    }
    int numPages;
    int maxEntriesPerPage;
    int pageSize;
    dataPage *page;
};

// ======================================== 
// the OPTICS algorithm
#define HANDLED -1
#define NOTHANDLED -2
// ======================================== 
// helper class to manage the seed list
class seedList {
public:
    seedList(int _size);
    ~seedList();
    void addRandomPoint(dataBase *db);
    void update       (int index, float64 r);
    void getMin       (int *index, float64 *r);
    int  getNumEntries(void)  { return(numEntries); }
private:
    void adjustHeap   (int index);
    int     numEntries;
    int     doneIndex;
    int     size;
    float64 * slReach;
    int *   slIndex;
    int *   slInvIndex;
};

// the optics algorithm
class optics {
public:
    optics               (dataBase *_db);
    ~optics              ();
    float64 updateAllReach (int index);
    void output          (int index, float64 coreLevel, float64 reach);
    //ptCluster **clusters;
    vector<ptCluster*> clusters;
    int nclusters;
private:
    dataBase *db;
    seedList *sl;
    float64 *ndist;
    char buffer[256];
};


