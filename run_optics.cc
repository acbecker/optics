#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <time.h>
#include "optics.H"
#include <math.h>
#include <vector>
#include <string>
using namespace std;

#define NS       100
#define SCRLEN   8096  /* those fits headers can be long! */
#define IDSIZE   64

/* defaults */
#define D_HTMDEPTH 25
#define D_TREECAP  50
#define D_TREEDIM  2
#define D_MINPTS   5
#define D_EPSILON  2

/* HTM */
extern uint64 cc_radec2ID(double ra, double dec, int depth);

//extern "C" {
//    int pix2wcs (struct WorldCoor *wcs,  /* World coordinate system structure */
//                 double xpix, 
//                 double ypix,            /* Image coordinates in pixels */
//                 double *xpos, 
//                 double *ypos            /* RA and Dec in degrees (returned) */
//                 );
//    struct WorldCoor *GetWCSFITS(char *filename, int verbose);
//    char *GetFITShead(char *filename, int verbose);
//};

/* global - these bastards caused so many problems... */
globalData GD;
int sortdim;

static void usage(void);
//extern struct WorldCoor *wcsxinit (double cra, double cdec, double secpix, double xrpix, double yrpix,
//                                   int nxpix, int nypix, double rotate, int equinox, double epoch, char *proj);


/* vectors : http://cplus.about.com/od/beginnerctutorial/l/aa050102a.htm */

int main (int argc, char ** argv) {
    char *filein, *val, line[SCRLEN];
    FILE *fptr;
    float64 x, y;
    uint64 htmid;
    long id;
    int i, j, nstars, nk=0, skiplines=0;
    string s, shead;
    
    vector<string> lines;
    vector<string> keywords;
    vector<float64> ras;
    vector<float64> decs;
    vector<uint64> htms;
    
    /* getopt */
    extern char *optarg;
    extern int optind;

    const float64 deg2rad = 0.01745329251994329576923690768488;
    
    /* HTM */
    int htmdepth = D_HTMDEPTH;
    
    /* tree */
    int dimension = D_TREEDIM; /* dimension of tree */
    int capacity  = D_TREECAP; /* number of points per index page */
    Data *d;
    
    /* optics */
    float64 epsilon = D_EPSILON; /* arcsec */ 
    int minpts     = D_MINPTS;
    ptCluster *pc;
    int dtype = 0;      /* euclidean distance measure */
    
    /* libwcs */
    //struct WorldCoor *WCS=NULL;
    float fakewcs = 0;  /* pixel scale of fake wcs projection */
    
    /* process command line options */
    while ((i = getopt(argc, argv, "3c:d:e:k:mn:s:w:")) != EOF) {
        switch (i) {
        case '3':
            dtype     = 1; /* use 3-d euclidean */
            dimension = 3;
            break;
        case 'c':
            capacity = atoi(optarg);
            break;
        case 'd':
            /* if htmdepth = 0, use ID from input file */
            htmdepth = atoi(optarg);
            break;
        case 'e':
            epsilon = atof(optarg);
            break;
        case 'k':
            for (val = strtok(optarg, ","); val; val=strtok(NULL, ","), nk++) {
                s = val;
                keywords.push_back(s.c_str());
                s.clear();
                fprintf(stderr, "Reading keyword %s\n", keywords.at(nk).c_str());
            }
            break;
        case 's':
            skiplines = atoi(optarg);
            break;
        case 'n':
            minpts = atoi(optarg);
            break;
        case 'w':
            fakewcs = atof(optarg);
            dimension = 2;
            dtype = 0;
            break;
        default:
            usage();
            exit(1);
        }
    }
    
    if (argc - optind < 1) {
        fprintf(stderr, "Please send infiles\n\n");
        usage();
        exit(1);
    }
    
    epsilon /= 3600;     /* in degrees */
    epsilon *= deg2rad;  /* in radians */
    
    if (dtype == 1) {
        /* angular separation = a
           v1 = (x1,y1,z2)
           v2 = (x2,y2,z2)
           v1.dot.v2 = cos(a)
           
           (x1 - x2)2 + (y1 - y2)2 + (z1 - z2)2 =
           v1.dot.v1 + v2.dot.v2 - 2(x1*x2 + y1*y2 + z1*z2) =
           2(1 - v1.dot.v2) = 2(1 - cos(a))
        */
        cout << "Clustering Radius : " << epsilon << " radians" << endl;
        epsilon  = sqrt(2. * (1. - cos(epsilon)));
        cout << "Clustering Radius : " << epsilon << " in sqrt(2*(1-cos(theta)))" << endl;
    } 
    else {
        cout << "Clustering Radius : " << epsilon << " radians" << endl;
    }
    
    /* set up GD */
    GD.dtype  = dtype;     /* what is the distance metric type */
    GD.dim    = dimension; /* number of dimensions to cluster on */
    GD.eps    = epsilon;   /* epsilon value for clustering algorithm */
    GD.minPts = minpts;    /* minPts value for clustering algorithm */
    
    /* set up index */
    /* THIS WAS IMPORTANT TO PUT THIS DOWN HERE! */
    Leanindex <Data> leanindex (dimension, 0);
    
    nstars = 0;
    
    time_t Tstart, Tend;
    double CPUtime;
    /* ********************** */
    /* Read in star positions */
    time (&Tstart);
    
    for (i = optind; i < argc; i++) {
        filein = argv[i];
        fptr   = fopen(filein, "r");
        
        shead  = filein;
        shead.append(" ");

        /* skip lines if requested */
        for (j = 0; j < skiplines; j++) 
            fgets(line, SCRLEN, fptr);
        
        while(fgets(line, SCRLEN, fptr) != NULL) {
            /* skip poorly formatted lines... */
            if (isalpha(line[0])) continue;
            if (htmdepth == 0) {
                if (sscanf(line, "%ld %lf %lf", &id, &x, &y) == 3) {
                    s.clear();
                    s.append(shead);  
                   
                    ras.push_back(x * deg2rad);  /* its actually ra  (DEGREES!) */
                    decs.push_back(y * deg2rad); /* its actually dec (DEGREES!) */
                    htms.push_back(id);

                    //fprintf(stdout, "#CAW1 %lf %.11e    %lf %.11e\n", x, ras.back(), y, decs.back());
                    
                    s.append(line);
                    lines.push_back(s.c_str());
                    nstars++;
                }
            }
            else {
                if (sscanf(line, "%lf %lf", &x, &y) == 2) {
                    s.clear();
                    s.append(shead);  
                    
                    ras.push_back(x * deg2rad);  /* its actually ra  (DEGREES!) */
                    decs.push_back(y * deg2rad); /* its actually dec (DEGREES!) */
                    htms.push_back(cc_radec2ID(x, y, htmdepth));
                    
                    s.append(line);
                    lines.push_back(s.c_str());
                    nstars++;
                }
            }
        }
        fclose(fptr);
        
    }
    
    
    time (&Tend);
    CPUtime= difftime(Tend, Tstart);
    printf("\n\n**** Performance\n");
    printf("     total cpu time used reading in objects : %.2lf sec.\n", CPUtime);
    /* Read in star positions */
    /* ********************** */
    
    /* ********************** */
    /* Fill up a tree         */
    time (&Tstart);
    for (i = nstars; i--; ) {
        d = new Data(GD.dim, i);
        if (dtype == 1) {
            d->point[0] = cos(decs.at(i)) * cos(ras.at(i));
            d->point[1] = cos(decs.at(i)) * sin(ras.at(i));
            d->point[2] = sin(decs.at(i));
            //if ( (i == 3) or (i == 5) )
            //fprintf(stdout, "#CAW2 %d %.11e %.11e %.11e %.11e %.11e\n", i, ras.at(i), decs.at(i), d->point[0], d->point[1], d->point[2]);
        }
        else {
            d->point[0] = ras.at(i);
            d->point[1] = decs.at(i);
        }
        leanindex.insert(d);
    }
    /* Somehow the data.id get screwed up here, damn... */
    leanindex.construct(capacity);
    
    /* if we want to dump the DB */
    //fptr = fopen("foo.dat", "w");
    //leanindex.out(fptr);
    //fclose(fptr);
    
    time (&Tend);
    CPUtime= difftime(Tend, Tstart);
    printf("\n\n**** Performance\n");
    printf("     total cpu time used filling and indexing tree : %.2lf sec.\n", CPUtime);
    /* Fill up a tree         */
    /* ********************** */
    
    /* ********************** */
    /* Now do some optics on it! */
    
    leanTreeDataBase db;
    if (db.load(leanindex) != 0)
        exit(2);
    
    printf("\n**** Executing OPTICS (based on lean tree)\n");
    time (&Tstart);
    
    optics *op = new optics(&db); 
    assert(op);
    time (&Tend);

    /*
    for (i = 0; i < 12; i++) {
        const point p0 = db.getPoint(i);

        for (j = i+1; j < 12; j++) {
            const point p1 = db.getPoint(j);

            cout << "#0 " << p0.getId() << " " << p0.getCoord(0) << " " << p0.getCoord(1) << " " << p0.getCoord(2) << " " << lines.at(p0.getId());
            cout << "#1 " << p1.getId() << " " << p1.getCoord(0) << " " << p1.getCoord(1) << " " << p1.getCoord(2) << " " << lines.at(p1.getId());
            cout << p0.dist(p1) << endl << endl;
        }
    }
    */

    CPUtime= difftime(Tend, Tstart);
    printf("\n\n**** Performance\n");
    printf("     total cpu time used running optics  : %.2lf sec.\n", CPUtime);
    
    cout.precision(6);
    
    cout << "# Found " << op->nclusters << " clusters" << endl;
    
    for (i = op->nclusters; i--; ) {
        pc = op->clusters.at(i);
        pc->center();
        if (pc->npts >= GD.minPts) {
            fprintf(stdout, "CLUSTER %d : %d %.6f %.6f\n", i, pc->npts, pc->ra/deg2rad, pc->dec/deg2rad);
            for (j = pc->npts; j--; ) {
                const point p = db.getPoint(pc->inds[j]);
                
                /* not really htmid, just running index */
                htmid = p.getId();
                cout << htms.at(htmid) << " " <<  lines.at(htmid);
            }
        }
    }
    
    delete op;
    
    /* Now do some optics on it! */
    /* ********************** */
    
    exit(0);
}


static void usage(void) {
    cout << "Usage: run_optics [options] infiles" << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << "          -3     : use 3-d euclidean distance metric" << endl;
    cout << "          -c int : number of points per index page in tree (" << D_TREECAP << ")" << endl;
    cout << "          -d int : HTM depth for RA,DEC storage (" << D_HTMDEPTH << ")" << endl;
    cout << "          -e flt : epsilon for clustering in arcsec (" << D_EPSILON << ")" << endl;
    cout << "          -h     : use Haversine distance metric (very slow)" << endl;
    cout << "          -k a,b : comma separated list of header keywords to keep, valid for .cmp" << endl;
    cout << "          -n int : minimum number of points in cluster (" << D_MINPTS << ")" << endl;
    cout << "          -s int : number of header lines to skip for each input file" << endl;
    //cout << "          -w flt : make fake TAN WCS projection to project points into, with this pixel scale" << endl;
}
