/*
 jit.boids3d.c -- Jack Truskowski and Grace Handler 2015-2016
 
 adapted from:
 12/15/2005 wesley smith
 boids3d 08/2005 a.sier / jasch adapted from boids by eric singer ÔøΩ 1995-2003 eric l. singer
 free for non-commercial use
 modified 070207 by sier to set boids pos, etc
 
 */

#include "jit.common.h"
#include <math.h>

// constants
#define			kAssistInlet	1
#define			kAssistOutlet	2
#define			kMaxLong		0xFFFFFFFF
#define			kMaxNeighbors	200

//Maximum number of flocks allowed in the simulation
#define         MAX_FLOCKS      6

// initial flight parameters
const int           kBoidMaxAge     = 1000; // default max age
const long			kNumBoids		= 0;	// number of boids for each flock
const long			kNumNeighbors	= 2;	// must be <= kMaxNeighbors
const double 		kMinSpeed		= 0.15;	// boids' minimum speed
const double		kMaxSpeed		= 0.25;	// boids' maximum speed
const double		kCenterWeight	= 0.25;	// flock centering
const double		kAttractWeight	= 0.300;// attraction point seeking
const double		kMatchWeight	= 0.100;// neighbors velocity matching
const double		kAvoidWeight	= 0.10;	// neighbors avoidance
const double		kWallsWeight	= 0.500;// wall avoidance [210] --  NOT USED
const double		kEdgeDist		= 0.5;	// vision distance to avoid wall edges [5]  --  NOT USED
const double		kDefaultSpeed	= 0.100;// Default boid speed
const double		kInertiaFactor	= 0.20;	// willingness to change speed & direction
const double		kAccelFactor	= 0.100;// neighbor avoidance accelerate or decelerate rate
const double        kNRadius        = 0.25; // neighborhood radius
const double		kFlyRectTop		= 1.0;	// fly rect boundaries
const double		kFlyRectLeft	= -1.0;
const double		kFlyRectBottom	= -1.0;
const double		kFlyRectRight	= 1.0;
const double		kFlyRectFront	= 1.0;
const double		kFlyRectBack	= -1.0;
const double        kFlyRectScalingFactor = 10;

//use defines instead of structs in jitter object
//because of the way attribute data types work

//defines for Point3d and Velocity
#define				x			0
#define				y			1
#define				z			2

//defines for FlyRect
#define				left		0
#define				right		1
#define				top			2
#define				bottom		3
#define				front		4
#define				back		5

/*
    Struct for an attract pt
 */
typedef struct Attractor {
    
    int id;
    
    double loc[3];
    
    struct Attractor *nextAttractor;
    
} Attractor, *AttractorPtr;

/*
    Struct for a boid object
        Boids are stored in LinkedLists by flock in the _jit_boids3d struct
 */
typedef struct Boid {
    
    int         flockID;
    int         age;
    
	double		oldPos[3];
	double		newPos[3];
	double		oldDir[3];
	double		newDir[3];
	double		speed;
	long		neighbor[kMaxNeighbors];
	double		neighborDistSqr[kMaxNeighbors];
    
    struct Boid *nextBoid;
    
} Boid, *BoidPtr;

/*
    Struct for the actual jitter object
 */
typedef struct _jit_boids3d
{
	t_object		ob;
	char			mode;
	long			number;
    long            numAttractors;
	long			neighbors;
	
	double			flyrect[6]; //dimensions of the simulation
	long			flyRectCount;
    
    int             allowNeighborsFromDiffFlock; //can boids find neighbors that aren't in their flock
	
    //Parameters for each flock
    int             boidCount[MAX_FLOCKS];
    int             flockID[MAX_FLOCKS]; //?
	double 			minspeed[MAX_FLOCKS];
	double			maxspeed[MAX_FLOCKS];
	double			center[MAX_FLOCKS];
	double			attract[MAX_FLOCKS];
	double			match[MAX_FLOCKS];
	double			avoid[MAX_FLOCKS];
	double			repel[MAX_FLOCKS];
	double			edgedist[MAX_FLOCKS];
	double			speed[MAX_FLOCKS];
	double			inertia[MAX_FLOCKS];
	double			accel[MAX_FLOCKS];
    double          neighborRadius[MAX_FLOCKS];
    double          age[MAX_FLOCKS];
	
    //XYZ for each flock
    
    double          tempCenterPt[3];
	long			centerPtCount;
	
    BoidPtr         flockLL[MAX_FLOCKS]; //array holding 6 linked lists of the flocks
    AttractorPtr    attractorLL; //linked list for attractors
	
	double 			d2r;
	double			r2d;
    
} t_jit_boids3d;


/*
    Methods for the jitter object
 */
void *_jit_boids3d_class;
t_jit_err jit_boids3d_init(void);
t_jit_boids3d *jit_boids3d_new(void);
void freeFlocks(t_jit_boids3d *flockPtr);
t_jit_err jit_boids3d_matrix_calc(t_jit_boids3d *flockPtr, void *inputs, void *outputs);

void jit_boids3d_calculate_ndim(t_jit_boids3d *flockPtr, long dimcount, long *dim, long planecount,
                                t_jit_matrix_info *out_minfo, char *bop);

//Attribute methods
/*
    These are the method that will get called when a message is received from the Max patch
        ie) [pak age 0. 0.] -> will call the jit_boids3d_age() method
 */
t_jit_err jit_boids3d_neighbors(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
t_jit_err jit_boids3d_minspeed(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
t_jit_err jit_boids3d_nradius(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
t_jit_err jit_boids3d_inertia(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
t_jit_err jit_boids3d_number(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
t_jit_err jit_boids3d_maxspeed(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
t_jit_err jit_boids3d_center(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
t_jit_err jit_boids3d_attract(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
t_jit_err jit_boids3d_match(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
t_jit_err jit_boids3d_avoid(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
t_jit_err jit_boids3d_repel(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
t_jit_err jit_boids3d_edgedist(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv); //NOT USED
t_jit_err jit_boids3d_speed(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
t_jit_err jit_boids3d_accel(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
t_jit_err jit_boids3d_age(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
t_jit_err jit_boids3d_attractpt(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
t_jit_err jit_boids3d_addattractor(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
t_jit_err jit_boids3d_deleteattractor(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);


//Initialization methods
void InitFlock(t_jit_boids3d *flockPtr);
BoidPtr InitLL(t_jit_boids3d *flockPtr, long numBoids, int flockID); //void Flock_donumBoids(t_jit_boids3d *flockPtr, long numBoids);
BoidPtr InitBoid(t_jit_boids3d *flockPtr);
AttractorPtr InitAttractor(t_jit_boids3d *flockPtr);

//Methods for running the simulation
void FlightStep(t_jit_boids3d *flockPtr);
float CalcFlockCenterAndNeighborVel(t_jit_boids3d *flockPtr, BoidPtr theBoid, double *matchNeighborVel, double *avoidNeighborVel);
void SeekPoint(t_jit_boids3d *flockPtr, BoidPtr theBoid, double *seekPt, double* seekDir);
void SeekAttractors(t_jit_boids3d *flockPtr, BoidPtr theBoid, double* seekDir);
void AvoidWalls(t_jit_boids3d *flockPtr, BoidPtr theBoid, double *wallVel);
char InFront(BoidPtr theBoid, BoidPtr neighbor);
void NormalizeVelocity(double *direction);
double RandomInt(double minRange, double maxRange);
double DistSqrToPt(double *firstPoint, double *secondPoint);
int CalcNumBoids(t_jit_boids3d *flockPtr);


/*
    Initializes the jitter object
 */
t_jit_err jit_boids3d_init(void)
{
    long attrflags=0;
	t_jit_object *attr,*mop,*o, *o2, *o3; //o, o2, and o3 are the 3 outlets. Mop stands for a matrix in jitter
	t_symbol *atsym;
	
	atsym = gensym("jit_attr_offset");
	
    //make a class and tell it the methods it will use to initialize and free
	_jit_boids3d_class = jit_class_new("jit_boids3d",(method)jit_boids3d_new,(method)freeFlocks,
                                       sizeof(t_jit_boids3d),0L);
    
	//add mop
	mop = jit_object_new(_jit_sym_jit_mop,0,3); //object will have 0 inlets and 3 outlets
	o = jit_object_method(mop,_jit_sym_getoutput,1); //first outlet
    o2 = jit_object_method(mop,_jit_sym_getoutput,2); //second outlet
    o3 = jit_object_method(mop,_jit_sym_getoutput,3); //second outlet
	jit_attr_setlong(o,_jit_sym_dimlink,0);
    jit_attr_setlong(o2,_jit_sym_dimlink,0);
    jit_attr_setlong(o3,_jit_sym_dimlink,0);

	
	jit_class_addadornment(_jit_boids3d_class,mop);
	//add methods
	jit_class_addmethod(_jit_boids3d_class, (method)jit_boids3d_matrix_calc, 		"matrix_calc", 		A_CANT, 0L);
	jit_class_addmethod(_jit_boids3d_class, (method)InitBoid, 				"init_boid", 			A_USURP_LOW, 0L);
    
	//add attributes
	attrflags = JIT_ATTR_GET_DEFER_LOW | JIT_ATTR_SET_USURP_LOW;
    
    //mode
	attr = jit_object_new(atsym,"mode",_jit_sym_char,attrflags,
                          (method)0L,(method)0L,calcoffset(t_jit_boids3d,mode));
	jit_class_addattr(_jit_boids3d_class,attr);
    
    //allow boids from diff flocks
	attr = jit_object_new(atsym,"diffFlock",_jit_sym_char,attrflags,
                          (method)0L,(method)0L,calcoffset(t_jit_boids3d,allowNeighborsFromDiffFlock));
	jit_class_addattr(_jit_boids3d_class,attr);
    
    //neighbor radius
	attr = jit_object_new(_jit_sym_jit_attr_offset_array,"nradius",_jit_sym_float64,2,attrflags,
                          (method)0L,(method)jit_boids3d_nradius, calcoffset(t_jit_boids3d,neighborRadius));
	jit_class_addattr(_jit_boids3d_class,attr);
    
	//number
	attr = jit_object_new(_jit_sym_jit_attr_offset_array,"number",_jit_sym_long, 6, attrflags,
                          (method)0L,(method)jit_boids3d_number,calcoffset(t_jit_boids3d,number));
	jit_class_addattr(_jit_boids3d_class,attr);
	
	//neighbors
	attr = jit_object_new(atsym,"neighbors",_jit_sym_long,attrflags,
                          (method)0L,(method)jit_boids3d_neighbors,calcoffset(t_jit_boids3d,neighbors));
	jit_class_addattr(_jit_boids3d_class,attr);
    
	//flyrect
	attr = jit_object_new(_jit_sym_jit_attr_offset_array,"flyrect",_jit_sym_float64,6,attrflags,
                          (method)0L,(method)0L,calcoffset(t_jit_boids3d,flyRectCount),calcoffset(t_jit_boids3d,flyrect));
	jit_class_addattr(_jit_boids3d_class,attr);
	
	//minspeed
	attr = jit_object_new(_jit_sym_jit_attr_offset_array,"minspeed",_jit_sym_float64,2,attrflags,
                          (method)0L,(method)jit_boids3d_minspeed,calcoffset(t_jit_boids3d,minspeed));
	jit_class_addattr(_jit_boids3d_class,attr);
	
	//maxspeed
	attr = jit_object_new(_jit_sym_jit_attr_offset_array,"maxspeed",_jit_sym_float64,2,attrflags,
                          (method)0L,(method)jit_boids3d_maxspeed,calcoffset(t_jit_boids3d,maxspeed));
	jit_class_addattr(_jit_boids3d_class,attr);
	
	//center
	attr = jit_object_new(_jit_sym_jit_attr_offset_array,"center",_jit_sym_float64,2,attrflags,
                          (method)0L,(method)jit_boids3d_center,calcoffset(t_jit_boids3d,center));
	jit_class_addattr(_jit_boids3d_class,attr);
	
	//attract
	attr = jit_object_new(_jit_sym_jit_attr_offset_array,"attract",_jit_sym_float64,2,attrflags,
                          (method)0L,(method)jit_boids3d_attract,calcoffset(t_jit_boids3d,attract));
	jit_class_addattr(_jit_boids3d_class,attr);
	
	//match
	attr = jit_object_new(_jit_sym_jit_attr_offset_array,"match",_jit_sym_float64,2,attrflags,
                          (method)0L,(method)jit_boids3d_match,calcoffset(t_jit_boids3d,match));
	jit_class_addattr(_jit_boids3d_class,attr);
	
	//avoid
	attr = jit_object_new(_jit_sym_jit_attr_offset_array,"avoid",_jit_sym_float64,2,attrflags,
                          (method)0L,(method)jit_boids3d_avoid,calcoffset(t_jit_boids3d,avoid));
	jit_class_addattr(_jit_boids3d_class,attr);
	
	//repel
	attr = jit_object_new(_jit_sym_jit_attr_offset_array,"repel",_jit_sym_float64,2,attrflags,
                          (method)0L,(method)jit_boids3d_repel,calcoffset(t_jit_boids3d,repel));
	jit_class_addattr(_jit_boids3d_class,attr);
	
	//edgedist
	attr = jit_object_new(_jit_sym_jit_attr_offset_array,"edgedist",_jit_sym_float64,2,attrflags,
                          (method)0L,(method)jit_boids3d_edgedist,calcoffset(t_jit_boids3d,edgedist));
	jit_class_addattr(_jit_boids3d_class,attr);
	
	//speed
	attr = jit_object_new(_jit_sym_jit_attr_offset_array,"speed",_jit_sym_float64,2,attrflags,
                          (method)0L,(method)jit_boids3d_speed,calcoffset(t_jit_boids3d,speed));
	jit_class_addattr(_jit_boids3d_class,attr);
	
	//inertia
	attr = jit_object_new(_jit_sym_jit_attr_offset_array,"inertia",_jit_sym_float64,2,attrflags,
                          (method)0L,(method)jit_boids3d_inertia,calcoffset(t_jit_boids3d,inertia));
	jit_class_addattr(_jit_boids3d_class,attr);
	
	//accel
	attr = jit_object_new(_jit_sym_jit_attr_offset_array,"accel",_jit_sym_float64,2,attrflags,
                          (method)0L,(method)jit_boids3d_accel,calcoffset(t_jit_boids3d,accel));
	jit_class_addattr(_jit_boids3d_class,attr);
    
    //attractpt
    attr = jit_object_new(_jit_sym_jit_attr_offset_array,"attractpt",_jit_sym_float64, 4, attrflags,
                          (method)0L,(method)jit_boids3d_attractpt,calcoffset(t_jit_boids3d,numAttractors));
    jit_class_addattr(_jit_boids3d_class,attr);
    
    //age
    attr = jit_object_new(_jit_sym_jit_attr_offset_array,"age",_jit_sym_float64,2,attrflags,
                          (method)0L,(method)jit_boids3d_age,calcoffset(t_jit_boids3d,speed));
    jit_class_addattr(_jit_boids3d_class,attr);
    
    //add attractor
    attr = jit_object_new(_jit_sym_jit_attr_offset_array,"addattractor",_jit_sym_long,2,attrflags,
                          (method)0L,(method)jit_boids3d_addattractor,calcoffset(t_jit_boids3d,numAttractors));
    jit_class_addattr(_jit_boids3d_class,attr);
    
    //add attractor
    attr = jit_object_new(_jit_sym_jit_attr_offset_array,"deleteattractor",_jit_sym_long,2,attrflags,
                          (method)0L,(method)jit_boids3d_deleteattractor,calcoffset(t_jit_boids3d,numAttractors));
    jit_class_addattr(_jit_boids3d_class,attr);
    
	
	jit_class_register(_jit_boids3d_class); //register the class with Max
    
	return JIT_ERR_NONE;
    
}

//boids attribute methods

/*
    The following methods are for moving/adding/deleting attractors
 */

/*
    Updates the position of an attractor
    Inputs: 
        argv+0/1/2 = xyz position to update to
        argv+3 = ID of the attractor
 */
t_jit_err jit_boids3d_attractpt(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
    int attractorID = (int)jit_atom_getfloat(argv+3);
    AttractorPtr iterator = flockPtr->attractorLL;
    while (iterator){
        if(attractorID == iterator->id){
            //this is the attractor we want to modify
            if(iterator){
                iterator->loc[0] = (double)jit_atom_getfloat(argv);
                iterator->loc[1] = (double)jit_atom_getfloat(argv+1);
                iterator->loc[2] = (double)jit_atom_getfloat(argv+2);
            }
            return JIT_ERR_NONE;
        }
        iterator = iterator->nextAttractor;
    }
    
    return JIT_ERR_NONE;
}


/*
    Adds an attractor at the origin
    Inputs: none
 */
t_jit_err jit_boids3d_addattractor(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
    //grab the ID of the new attractor
    int newID = (int)jit_atom_getlong(argv);
    
    //initialize an attractor
    AttractorPtr iterator = flockPtr->attractorLL;
    AttractorPtr newAttractor = InitAttractor(flockPtr);
    flockPtr->numAttractors++;
    
    //no attractors exist
    if(!iterator){
        newAttractor->id = newID;
        flockPtr->attractorLL = newAttractor;
        flockPtr->attractorLL->nextAttractor = NULL;
        return JIT_ERR_NONE;
    }
    
    //at least one attractor already exists
    //check if there already exists an attractor with newID
    int idAlreadyExists = 0;
    int maxID = 0;
    while(iterator){
        if(iterator->id > maxID){
            maxID = iterator->id;
        }
        if(iterator->id == newID){
            idAlreadyExists = 1;
            break;
        }
        iterator = iterator->nextAttractor;
    }
    iterator = flockPtr->attractorLL;
    
    //add to front
    if(idAlreadyExists == 0){
        newAttractor->id = newID;
    }else{
        newAttractor->id = maxID+1;
    }
    newAttractor->nextAttractor = iterator;
    flockPtr->attractorLL = newAttractor;
    return JIT_ERR_NONE;
}


/*
    Deletes a specific attractor
    Inputs: 
        argv+0: ID of the attractor to be deleted
 */
t_jit_err jit_boids3d_deleteattractor(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
    int attractorID = (int)jit_atom_getlong(argv);
    AttractorPtr iterator = flockPtr->attractorLL;
    AttractorPtr prev = iterator;
    
    //iterate thru the LL of attractors to find the right one
    while(iterator){
        if(iterator->id == attractorID){
            //this is the attractor to delete
            if(prev == iterator){ //first attractor
                iterator = iterator->nextAttractor;
                free(prev);
                flockPtr->attractorLL = iterator;
                
            }else{
                
                //middle or end
                prev->nextAttractor = iterator->nextAttractor;
                free(iterator);
            }
            
            //update numAttractors and mark the LL NULL if necessary
            flockPtr->numAttractors--;
            if(flockPtr->numAttractors <= 0){
                flockPtr->attractorLL = NULL;
            }
            return JIT_ERR_NONE;
            
        }
        prev = iterator;
        iterator = iterator->nextAttractor;
    }
    
    //couldn't find this attractor
    return JIT_ERR_NONE;
}


/*
    Following methods update flock-specific attributes
 */

 //---NOT CURRENTLY USED---
t_jit_err jit_boids3d_neighbors(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
    flockPtr->neighbors = (double)MIN(jit_atom_getfloat(argv), kMaxNeighbors);
    return JIT_ERR_NONE;
}
t_jit_err jit_boids3d_nradius(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
    int flockID = (int)jit_atom_getfloat(argv+1);
	flockPtr->neighborRadius[flockID] = (double)MAX(jit_atom_getfloat(argv), 0.0);
    return JIT_ERR_NONE;
}

t_jit_err jit_boids3d_minspeed(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
    int flockID = (int)jit_atom_getfloat(argv+1);
	flockPtr->minspeed[flockID] = (double)MAX(jit_atom_getfloat(argv), 0.000001);
    return JIT_ERR_NONE;
}

t_jit_err jit_boids3d_maxspeed(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
    int flockID = (int)jit_atom_getfloat(argv+1);
	flockPtr->maxspeed[flockID] = (double)MAX(jit_atom_getfloat(argv), 0.000001);
    return JIT_ERR_NONE;
}

t_jit_err jit_boids3d_center(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
    int flockID = (int)jit_atom_getfloat(argv+1);
	flockPtr->center[flockID] = (double)MAX(jit_atom_getfloat(argv), 0.000001);
    return JIT_ERR_NONE;
}

t_jit_err jit_boids3d_attract(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
    int flockID = (int)jit_atom_getfloat(argv+1);
	flockPtr->attract[flockID] = (double)MAX(jit_atom_getfloat(argv), 0.000001);
    return JIT_ERR_NONE;
}

t_jit_err jit_boids3d_match(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
    int flockID = (int)jit_atom_getfloat(argv+1);
	flockPtr->match[flockID] = (double)MAX(jit_atom_getfloat(argv), 0.000001);
    return JIT_ERR_NONE;
}

t_jit_err jit_boids3d_avoid(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
    int flockID = (int)jit_atom_getfloat(argv+1);
	flockPtr->avoid[flockID] = (double)MAX(jit_atom_getfloat(argv), 0.000001);
    return JIT_ERR_NONE;
}

t_jit_err jit_boids3d_repel(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
    int flockID = (int)jit_atom_getfloat(argv+1);
	flockPtr->repel[flockID] = (double)MAX(jit_atom_getfloat(argv), 0.000001);
    return JIT_ERR_NONE;
}

t_jit_err jit_boids3d_edgedist(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
    int flockID = (int)jit_atom_getfloat(argv+1);
	flockPtr->edgedist[flockID] = (double)MAX(jit_atom_getfloat(argv), 0.000001);
    return JIT_ERR_NONE;
}

t_jit_err jit_boids3d_speed(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
    int flockID = (int)jit_atom_getfloat(argv+1);
	flockPtr->speed[flockID] = (double)MAX(jit_atom_getfloat(argv), 0.000001);
    return JIT_ERR_NONE;
}

t_jit_err jit_boids3d_inertia(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
	double val = (double)jit_atom_getfloat(argv);
    int flockID =(int)jit_atom_getfloat(argv+1);
	
	if(val == 0.0)
		flockPtr->inertia[flockID] = 0.000001;
	else
		flockPtr->inertia[flockID] = val;
    return JIT_ERR_NONE;
}

t_jit_err jit_boids3d_accel(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
    int flockID = (int)jit_atom_getfloat(argv+1);
	flockPtr->accel[flockID] = (double)MAX(jit_atom_getfloat(argv), 0.000001);
    return JIT_ERR_NONE;
}

t_jit_err jit_boids3d_age(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
    int flockID = (int)jit_atom_getfloat(argv+1);
    flockPtr->age[flockID] = (double)jit_atom_getfloat(argv);
    return JIT_ERR_NONE;
}


/*
 This function is called when the number of boids in the maxpatch is changed. Adds or deletes the specified number of boids from the flock LL
 */
t_jit_err jit_boids3d_number(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
    int boidChanges[6]; //number of boids being deleted from each flock
    
    int newNumBoids; //new total number of boids across all flocks
    newNumBoids = jit_atom_getlong(argv+6);
    
    for (int i=0; i<6; i++){
        boidChanges[i] = (int)jit_atom_getlong(argv+i); //new total boids for the ith flock
    }
    
    //make sure the number of boids in at least one flock is being changed
    int changed = 0;
    for(int i=0; i<6; i++){
        if (boidChanges[i] != 0){
            changed = 1;
            break;
        }
    }
    if(changed == 0){
        return NULL;
    }
    
    //iterate thru flocks and update boids
    for (int i=0; i<MAX_FLOCKS; i++){
        
        if(boidChanges[i] == 0) continue; //no changes in this flock
        
        else if(boidChanges[i] < 0){ //we're deleting boids
            BoidPtr iterator = flockPtr->flockLL[i];
            BoidPtr toBeDeleted = iterator;
            while (boidChanges[i] < 0){
                if(!toBeDeleted){
                    break;
                }
                iterator = iterator->nextBoid;
                free(toBeDeleted);
                toBeDeleted = iterator;
                flockPtr->flockLL[i] = iterator;
                flockPtr->boidCount[i]--; //update the number of boids in flock
                boidChanges[i]++;
            }
        }else{ //we're adding boids
            for (int j=0; j<boidChanges[i]; j++){
                
                //initialize a new boid and add it to the front of the LL
                BoidPtr newBoid = InitBoid(flockPtr);
                if(!newBoid){
                    exit(1);
                }
                newBoid->nextBoid = flockPtr->flockLL[i];
                flockPtr->flockLL[i] = newBoid;
                newBoid->flockID = i;
                flockPtr->boidCount[i]++; //update the number of boids in flock
            }
        }
    }
    
    return 0;
}


/*
 Prepares the output matrix and sends it back to the max patch
 */
t_jit_err jit_boids3d_matrix_calc(t_jit_boids3d *flockPtr, void *inputs, void *outputs)
{
    
	t_jit_err err=JIT_ERR_NONE;
	long out_savelock, out2_savelock, out3_savelock; //if there is a problem, saves and locks the output matricies
	t_jit_matrix_info out_minfo, out2_minfo, out3_minfo;
	char *out_bp, *out2_bp, *out3_bp;
	long i,dimcount,planecount,dim[JIT_MATRIX_MAX_DIMCOUNT]; //dimensions and planes for the first output matrix
	void *out_matrix, *out2_matrix, *out3_matrix;
	
	out_matrix 	= jit_object_method(outputs,_jit_sym_getindex,0);
    out2_matrix 	= jit_object_method(outputs,_jit_sym_getindex,1);
    out3_matrix     = jit_object_method(outputs, _jit_sym_getindex, 2);
    
	if (flockPtr&&out_matrix&&out2_matrix&&out3_matrix) {
		out_savelock = (long) jit_object_method(out_matrix,_jit_sym_lock,1);
        out2_savelock = (long) jit_object_method(out2_matrix,_jit_sym_lock,1);
        out3_savelock = (long) jit_object_method(out3_matrix, _jit_sym_lock,1);
		
		jit_object_method(out_matrix,_jit_sym_getinfo,&out_minfo); //assign the out_infos to their cooresponding out matrix
        jit_object_method(out2_matrix,_jit_sym_getinfo,&out2_minfo);
        jit_object_method(out3_matrix,_jit_sym_getinfo, &out3_minfo);
		
        int numBoids = CalcNumBoids(flockPtr);
        
        //dimensions of the output matrix (number of boids x 1)
        out_minfo.dim[0] = numBoids;
		out_minfo.dim[1] = 1;
		out_minfo.type = _jit_sym_float32; //outputting floating point numbers
        
        //(number of flocks x 1)
        out2_minfo.dim[0] = MAX_FLOCKS;
		out2_minfo.dim[1] = 1;
		out2_minfo.type = _jit_sym_float32; //outputting floating point numbers
        out2_minfo.planecount = 1;
        
        //dimensions of attractor output matrix (number of attractors x 1)
        out3_minfo.dim[0] = flockPtr->numAttractors;
        out3_minfo.dim[1] = 1;
        out3_minfo.type = _jit_sym_float32; //outputting floating point numbers
        out3_minfo.planecount = 4; //xyz, id
        
        //output the correct mode
		switch(flockPtr->mode) { // newpos
			case 0:
				out_minfo.planecount = 4;
                break;
			case 1: //newpos + oldpos
				out_minfo.planecount = 7;
                break;
			case 2:	//newpos +  oldpos + speed-azimuth-elevation
				out_minfo.planecount = 10;
                break;
		}
		
        //for some reason, 2 of these in a row are required
		jit_object_method(out_matrix,_jit_sym_setinfo,&out_minfo);
		jit_object_method(out_matrix,_jit_sym_getinfo,&out_minfo);
        
        jit_object_method(out2_matrix,_jit_sym_setinfo,&out2_minfo);
		jit_object_method(out2_matrix,_jit_sym_getinfo,&out2_minfo);
        
        jit_object_method(out3_matrix,_jit_sym_setinfo,&out3_minfo);
        jit_object_method(out3_matrix,_jit_sym_getinfo,&out3_minfo);
		
		jit_object_method(out_matrix,_jit_sym_getdata,&out_bp);
        jit_object_method(out2_matrix,_jit_sym_getdata,&out2_bp);
        jit_object_method(out3_matrix,_jit_sym_getdata,&out3_bp);
        
        //something went wrong, handle the error
        if (!out_bp || !out2_bp || !out3_bp) {
            err=JIT_ERR_INVALID_OUTPUT;
            goto out;
        }
        
        //populate the second outlet matrix with data
        float *out2_data = (float*)out2_bp;
        for(int i=0; i<MAX_FLOCKS; i++){
            out2_data[0] = flockPtr->boidCount[i];
            out2_data+=1;
        }
        
        //populate the 3rd outlet with data
        float *out3_data = (float*)out3_bp;
        AttractorPtr iterator = flockPtr->attractorLL;
        while(iterator){
            out3_data[0] = iterator->loc[0];
            out3_data[1] = iterator->loc[1];
            out3_data[2] = iterator->loc[2];
            out3_data[3] = iterator->id;
        
            out3_data += 4; //planecount
            
            iterator=iterator->nextAttractor;
        }
        
        //get dimensions/planecount
        dimcount   = out_minfo.dimcount;
        planecount = out_minfo.planecount;
        
        for (i=0;i<dimcount;i++) {
            dim[i] = out_minfo.dim[i];
		}
        
        //populate the first outlet matrix with data
		jit_boids3d_calculate_ndim(flockPtr, dimcount, dim, planecount, &out_minfo, out_bp);
        
	} else {
		return JIT_ERR_INVALID_PTR;
	}
	
out: //output the matrix
	jit_object_method(out_matrix,gensym("lock"),out_savelock);
	return err;
}


/*
    Does a flight step and populates the first outlet matrix with the data (boids x,y,z etc)
 */
void jit_boids3d_calculate_ndim(t_jit_boids3d *flockPtr, long dimcount, long *dim, long planecount,
                                t_jit_matrix_info *out_minfo, char *bop)
{
	float *fop;
	double 	tempNew_x, tempNew_y, tempNew_z;
	double 	tempOld_x, tempOld_y, tempOld_z;
	double	delta_x, delta_y, delta_z, azi, ele, speed;
	
    //do a step in the simulation
	FlightStep(flockPtr);
	
	fop = (float *)bop; //contains the data
	
    //pick the correct mode (to get the appropriate number of planes)
    switch(flockPtr->mode) { // newpos
        case 0:
            for (int i=0; i<MAX_FLOCKS; i++){
                BoidPtr iterator = flockPtr->flockLL[i];
                
                while (iterator){ //iterate thru the boids in the flock and add their info to the matrix
                    fop[0] = iterator->newPos[x];
                    fop[1] = iterator->newPos[y];
                    fop[2] = iterator->newPos[z];
                    fop[3] = iterator->flockID;
                    
                    fop += planecount;
                    
                    iterator = iterator->nextBoid;
                }
            }
            break;
        case 1:
            for (int i=0; i<MAX_FLOCKS; i++){
                BoidPtr iterator = flockPtr->flockLL[i];
                
                while (iterator){ //iterate thru the boids in the flock and add their info to the matrix
                    fop[0] = iterator->newPos[x];
                    fop[1] = iterator->newPos[y];
                    fop[2] = iterator->newPos[z];
                    fop[3] = iterator->flockID;
                    fop[4] = iterator->oldPos[x];
                    fop[5] = iterator->oldPos[y];
                    fop[6] = iterator->oldPos[z];
                    
                    fop += planecount;
                    
                    iterator = iterator->nextBoid;
                }
            }
            break;
        case 2:
            for (int i=0; i<MAX_FLOCKS; i++){
                BoidPtr iterator = flockPtr->flockLL[i];
                
                while (iterator){ //iterate thru the boids in the flock and add their info to the matrix
                    tempNew_x = iterator->newPos[x];
                    tempNew_y = iterator->newPos[y];
                    tempNew_z = iterator->newPos[z];
                    tempOld_x = iterator->oldPos[x];
                    tempOld_y = iterator->oldPos[y];
                    tempOld_z = iterator->oldPos[z];
                    
                    delta_x = tempNew_x - tempOld_x;
                    delta_y = tempNew_y - tempOld_y;
                    delta_z = tempNew_z - tempOld_z;
                    azi = jit_math_atan2(delta_z, delta_x) * flockPtr->r2d;
                    ele = jit_math_atan2(delta_y, delta_x) * flockPtr->r2d;
                    speed = jit_math_sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);
                    
                    fop[0] = tempNew_x;
                    fop[1] = tempNew_y;
                    fop[2] = tempNew_z;
                    fop[3] = iterator->flockID;
                    fop[4] = tempOld_x;
                    fop[5] = tempOld_y;
                    fop[6] = tempOld_z;
                    fop[7] = speed;
                    fop[8] = azi;
                    fop[9] = ele;
                    
                    fop += planecount;
                    
                    iterator = iterator->nextBoid; //move to next boid
                }
            }
            break;
    }
}

/*
    MARK: ----------------------------Below methods do not relate to communication with the Max Patch, they are only used within the external---------------------------------
 */


/*
    Calculates each boid's new velocity and updates position
 */
void FlightStep(t_jit_boids3d *flockPtr)
{
    //All velocities start out at 0, otherwise weird boid velocities may happen
    double			goCenterVel[3] = {0,0,0};
    double			goAttractVel[3] = {0,0,0};
    double			matchNeighborVel[3] = {0,0,0};
    double			avoidNeighborVel[3] = {0,0,0};
    const double	zeroVel[3]	= {0.0, 0.0, 0.0};
    
    //get every boid from every flock
    for (int i=0; i<MAX_FLOCKS; i++){
        BoidPtr iterator = flockPtr->flockLL[i];
        BoidPtr prevBoid = NULL;
        
        while(iterator){ //grab every boid from this flock
            
            //update age and check if it's this boid's time to die
            iterator->age++;
            if(iterator->age > flockPtr->age[iterator->flockID] && flockPtr->age[iterator->flockID] != -1){
                
                BoidPtr deletor = iterator;
                
                //delete the boid and continue
                if(!prevBoid){ //this boid is at the head of the linked list
                    deletor = iterator;
                    iterator = iterator->nextBoid;
                    free(deletor);
                    flockPtr->flockLL[i] = iterator;
                    
                }else{ //this boid is somewhere in the middle of the LL
                    iterator = iterator->nextBoid;
                    free(deletor);
                    prevBoid->nextBoid = iterator;
                }
                
                //update boid pointers and count and move to next boid
                flockPtr->boidCount[i]--;
                continue;
                
            }
            
            //save position and velocity
            
            iterator->oldPos[x] = iterator->newPos[x];
            iterator->oldPos[y] = iterator->newPos[y];
            iterator->oldPos[z] = iterator->newPos[z];
            
            iterator->oldDir[x] = iterator->newDir[x];
            iterator->oldDir[y] = iterator->newDir[y];
            iterator->oldDir[z] = iterator->newDir[z];
            
            //calculate velocity updates
            int flockID = iterator->flockID;
            
            float avoidNeighborSpeed = CalcFlockCenterAndNeighborVel(flockPtr, iterator, matchNeighborVel,  avoidNeighborVel);
            
            //update velocity to include centering and attracting instincts
            SeekPoint(flockPtr, iterator, flockPtr->tempCenterPt, goCenterVel);
            
            //Seek the attractors
            SeekAttractors(flockPtr, iterator, goAttractVel);
            
            // compute resultant velocity using weights and inertia
            iterator->newDir[x] = flockPtr->inertia[flockID] * (iterator->oldDir[x]) +
            (flockPtr->center[flockID] * goCenterVel[x] +
             flockPtr->attract[flockID] * goAttractVel[x] +
             flockPtr->match[flockID] * matchNeighborVel[x] +
             flockPtr->avoid[flockID] * avoidNeighborVel[x]) / flockPtr->inertia[flockID];
            iterator->newDir[y] = flockPtr->inertia[flockID] * (iterator->oldDir[y]) +
            (flockPtr->center[flockID] * goCenterVel[y] +
             flockPtr->attract[flockID] * goAttractVel[y] +
             flockPtr->match[flockID] * matchNeighborVel[y] +
             flockPtr->avoid[flockID] * avoidNeighborVel[y]) / flockPtr->inertia[flockID];
            iterator->newDir[z] = flockPtr->inertia[flockID] * (iterator->oldDir[z]) +
            (flockPtr->center[flockID] * goCenterVel[z] +
             flockPtr->attract[flockID] * goAttractVel[z] +
             flockPtr->match[flockID] * matchNeighborVel[z] +
             flockPtr->avoid[flockID] * avoidNeighborVel[z]) / flockPtr->inertia[flockID];
            NormalizeVelocity(iterator->newDir);	// normalize velocity so its length is unified
            
            // set to avoidNeighborSpeed bounded by minspeed and maxspeed
            if ((avoidNeighborSpeed >= flockPtr->minspeed[flockID]) &&
                (avoidNeighborSpeed <= flockPtr->maxspeed[flockID]))
                iterator->speed = avoidNeighborSpeed;
            else if (avoidNeighborSpeed > flockPtr->maxspeed[flockID])
                iterator->speed = flockPtr->maxspeed[flockID];
            else
                iterator->speed = flockPtr->minspeed[flockID];
            
            
            //bounce back from walls if the boid is beyond the limit of the flyrect
            AvoidWalls(flockPtr, iterator, iterator->newDir);
            
            // calculate new position, applying speed
            iterator->newPos[x] += iterator->newDir[x] * iterator->speed * (flockPtr->speed[flockID] / 100.0);
            iterator->newPos[y] += iterator->newDir[y] * iterator->speed * (flockPtr->speed[flockID] / 100.0);
            iterator->newPos[z] += iterator->newDir[z] * iterator->speed * (flockPtr->speed[flockID] / 100.0);
            
            //move to next boid
            prevBoid = iterator;
            iterator = iterator->nextBoid;
            
        }
        
    }
}

/*
 Calculates the center of a flock, saves it in flockPtr->centerPt. Also computes avoid and matching of neighbor velocities
    >This is a method that is a combo of what was previously FindFlockCenter() and MatchAndAvoidNeighbors(), because both involve calculating the neighbors of a boid
 
 TODO: possible optimization = if we're only allowing boids from same flock, only calculate the center once for each flock
 */
float CalcFlockCenterAndNeighborVel(t_jit_boids3d *flockPtr, BoidPtr theBoid, double *matchNeighborVel, double *avoidNeighborVel)
{
    //Variables for centering
    int flockID = theBoid->flockID;
    double totalH = 0, totalV = 0, totalD = 0;
    
    //Matching
    matchNeighborVel[x] = 0;
    matchNeighborVel[y] = 0;
    matchNeighborVel[z] = 0;
    
    //Variables for avoidance
    double	avoidSpeed = theBoid->speed;
    int neighborsCount = 0; //counter to keep track of how many neighbors we've found
    
    float PreferredDistance = 0.1;
    
    for(int i=0; i<MAX_FLOCKS; i++){ //grab every boid
        
        BoidPtr iterator = flockPtr->flockLL[i];
        
        while(iterator){
            
            double dist = sqrt(DistSqrToPt(theBoid->oldPos, iterator->oldPos));
            
            if(dist < flockPtr->neighborRadius[flockID] && dist > 0.0 && neighborsCount < kMaxNeighbors){ //check if this boid is close enough to be a neighbor
                
                //check to ensure this boid is allowed / in same flock
                if (flockPtr->allowNeighborsFromDiffFlock == 0 && iterator->flockID != flockID){
                    continue;
                }
                
                //this boid is a neighbor
                
                //centering
                neighborsCount++;
                totalH += theBoid->oldPos[x];
                totalV += theBoid->oldPos[y];
                totalD += theBoid->oldPos[z];
                
                //matching
                matchNeighborVel[x] += iterator->oldDir[x];
                matchNeighborVel[y] += iterator->oldDir[y];
                matchNeighborVel[z] += iterator->oldDir[z];
                
                //
                
                if (InFront((theBoid), iterator)) {	// adjust speed
                    avoidSpeed /= (flockPtr->accel[flockID] / 100.0);
                }
                else {
                    avoidSpeed *= (flockPtr->accel[flockID] / 100.0);
                }
    
                neighborsCount++;
            }
            
            
            iterator = iterator->nextBoid; //move to next boid
        }
    }
    
    //normalize the matching velocity
    NormalizeVelocity(matchNeighborVel);
    
    //update the center point as an average of theBoid's neighbors
    if(neighborsCount > 0){ //get the average position of all boids in the flock
        flockPtr->tempCenterPt[x] = (double)	(totalH / neighborsCount);
        flockPtr->tempCenterPt[y] = (double)	(totalV / neighborsCount);
        flockPtr->tempCenterPt[z] = (double)	(totalD / neighborsCount);
    }else{ //only boid in flock, its position is the center point
        flockPtr->tempCenterPt[x] = theBoid->oldPos[x];
        flockPtr->tempCenterPt[y] = theBoid->oldPos[y];
        flockPtr->tempCenterPt[z] = theBoid->oldPos[z];
    }
    
    return avoidSpeed;
}


/*
    Calculates a new seek direction to a point
 */
void SeekPoint(t_jit_boids3d *flockPtr, BoidPtr theBoid, double *seekPt, double* seekDir)
{
	seekDir[x] = seekPt[x] - theBoid->oldPos[x];
	seekDir[y] = seekPt[y] - theBoid->oldPos[y];
	seekDir[z] = seekPt[z] - theBoid->oldPos[z];
	NormalizeVelocity(seekDir);
}


/*
    Calculates a new seek direction for all attractors
        >TODO: not sure if this works
 */
void SeekAttractors(t_jit_boids3d *flockPtr, BoidPtr theBoid, double* seekDir)
{
    double closestDir[3] = {0.0, 0.0, 0.0};
    double closestDist = DBL_MAX;
    AttractorPtr iterator = flockPtr->attractorLL;
    
    //iterate thru and sum up the direction to all attractors
    while(iterator){
        
        double dist = DistSqrToPt(iterator->loc, theBoid->oldPos);
        if(dist < closestDist){
            closestDist = dist;
            closestDir[0] = iterator->loc[0];
            closestDir[1] = iterator->loc[1];
            closestDir[2] = iterator->loc[2];
        }
        
        iterator = iterator->nextAttractor;
    }
    
    //update seekDir
    if(closestDist != DBL_MAX){
        seekDir[0] = closestDir[0] - theBoid->oldPos[0];
        seekDir[1] = closestDir[1] - theBoid->oldPos[1];
        seekDir[2] = closestDir[2] - theBoid->oldPos[2];
    }
    
    NormalizeVelocity(seekDir);
}


/*
    Bounces a boid back from the edges of the flyrect if necessary
 */
void AvoidWalls(t_jit_boids3d *flockPtr, BoidPtr theBoid, double *wallVel)
{
	double		testPoint[3];
    
    //	wallVel[x] = 0.0;
    //	wallVel[y] = 0.0;
    //	wallVel[z] = 0.0;
    
	/* calculate test point in front of the nose of the boid */
	/* distance depends on the boid's speed and the avoid edge constant */
	testPoint[x] = theBoid->oldPos[x] + theBoid->newDir[x] * (theBoid->speed * (flockPtr->speed[theBoid->flockID] / 100.0));// * flockPtr->edgedist[flockPtr->boid[theBoid].flockID];
    testPoint[y] = theBoid->oldPos[y] + theBoid->newDir[y] * (theBoid->speed * (flockPtr->speed[theBoid->flockID] / 100.0));// * flockPtr->edgedist[flockPtr->boid[theBoid].flockID];
    testPoint[z] = theBoid->oldPos[z] + theBoid->newDir[z] * (theBoid->speed * (flockPtr->speed[theBoid->flockID] / 100.0));// * flockPtr->edgedist[flockPtr->boid[theBoid].flockID];
    
    
	/* if test point is out of the left (right) side of flockPtr->flyrect, */
	/* return a positive (negative) horizontal velocity component */
	if (testPoint[x] < flockPtr->flyrect[left]*kFlyRectScalingFactor)
		wallVel[x] = ABS(wallVel[x]);
	else if (testPoint[x] > flockPtr->flyrect[right]*kFlyRectScalingFactor)
		wallVel[x] = - ABS(wallVel[x]);
    
	/* same with top and bottom */
	if (testPoint[y] > flockPtr->flyrect[top]*kFlyRectScalingFactor)
		wallVel[y] = - ABS(wallVel[y]);
	else if (testPoint[y] < flockPtr->flyrect[bottom]*kFlyRectScalingFactor)
		wallVel[y] = ABS(wallVel[y]);
    
	/* same with front and back*/
	if (testPoint[z] > flockPtr->flyrect[front]*kFlyRectScalingFactor)
		wallVel[z] = - ABS(wallVel[z]);
	else if (testPoint[z] < flockPtr->flyrect[back]*kFlyRectScalingFactor)
		wallVel[z] = ABS(wallVel[z]);
}


/*
    Determines if a neighbor boid is in front of another boid
 */
char InFront(BoidPtr theBoid, BoidPtr neighbor)
{
	float	grad, intercept;
	char result;
	
    /* we do this on 2 planes, xy, yz. if one returns false then we know its behind. a.sier/jasch 08/2005
     
     Find the gradient and y-intercept of a line passing through theBoid's oldPos
     perpendicular to its direction of motion.  Another boid is in front of theBoid
     if it is to the right or left of this linedepending on whether theBoid is moving
     right or left.  However, if theBoid is travelling vertically then just compare
     their vertical coordinates.
     
     */
	// xy plane
	
	// if theBoid is not travelling vertically...
	if (theBoid->oldDir[x] != 0) {
		// calculate gradient of a line _perpendicular_ to its direction (hence the minus)
		grad = -theBoid->oldDir[y] / theBoid->oldDir[x];
		
		// calculate where this line hits the y axis (from y = mx + c)
		intercept = theBoid->oldPos[y] - (grad * theBoid->oldPos[x]);
        
		/* compare the horizontal position of the neighbor boid with */
		/* the point on the line that has its vertical coordinate */
		if (neighbor->oldPos[x] >= ((neighbor->oldPos[y] - intercept) / grad)) {
			/* return true if the first boid's horizontal movement is +ve */
			result = (theBoid->oldDir[x] > 0);
            
			if (result==0) return 0;
			else goto next;
			
		} else {
			/* return true if the first boid's horizontal movement is +ve */
			result = (theBoid->oldDir[x] < 0);
			if (result==0) return 0;
			else goto next;
		}
	}
	/* else theBoid is travelling vertically, so just compare vertical coordinates */
	else if (theBoid->oldDir[y] > 0) {
		result = (neighbor->oldPos[y] > theBoid->oldPos[y]);
		if (result==0){
			return 0;
		}else{
			goto next;
		}
	}else{
		result = (neighbor->oldPos[y] < theBoid->oldPos[y]);
		if (result==0){
			return 0;
		} else {
			goto next;
		}
	}
next:
    
	// yz plane
	
	// if theBoid is not travelling vertically...
	if (theBoid->oldDir[y] != 0) {
		// calculate gradient of a line _perpendicular_ to its direction (hence the minus)
		grad = -theBoid->oldDir[z] / theBoid->oldDir[y];
		
		// calculate where this line hits the y axis (from y = mx + c)
		intercept = theBoid->oldPos[z] - (grad * theBoid->oldPos[y]);
        
		// compare the horizontal position of the neighbor boid with
		// the point on the line that has its vertical coordinate
		if (neighbor->oldPos[y] >= ((neighbor->oldPos[z] - intercept) / grad)) {
			// return true if the first boid's horizontal movement is +ve
			result = (theBoid->oldDir[y] > 0);
			if (result==0){
				return 0;
			}else{
				goto next2;
			}
		} else {
			// return true if the first boid's horizontal movement is +ve
			result = (theBoid->oldDir[y] < 0);
			if (result==0){
				return 0;
			}else{
				goto next2;
			}
		}
	}
	// else theBoid is travelling vertically, so just compare vertical coordinates
	else if (theBoid->oldDir[z] > 0) {
		result = (neighbor->oldPos[z] > theBoid->oldPos[z]);
		if (result==0){
			return 0;
		}else{
			goto next2;
		}
	}else{
		result = (neighbor->oldPos[z] < theBoid->oldPos[z]);
		if (result==0){
			return 0;
		}else{
			goto next2;
		}
	}
next2:
	return 1;
}

void NormalizeVelocity(double *direction)
{
	float	hypot;
	
	hypot = jit_math_sqrt(direction[x] * direction[x] + direction[y] * direction[y] + direction[z] * direction[z] );
    
	if (hypot != 0.0) {
		direction[x] = direction[x] / hypot;
		direction[y] = direction[y] / hypot;
		direction[z] = direction[z] / hypot;
	}
}

//Returns a random integer in a range
double RandomInt(double minRange, double maxRange)
{
	double	t, result;
	
	t = (double)(jit_rand() & 0x0000FFFF)/(double)(0x0000FFFF);
    
	result = (t * (maxRange - minRange)) + minRange;
	return(result);
}

double DistSqrToPt(double *firstPoint, double *secondPoint)
{
	double	a, b,c;
	a = firstPoint[x] - secondPoint[x];
	b = firstPoint[y] - secondPoint[y];
	c = firstPoint[z] - secondPoint[z];
	return(a * a + b * b + c * c);
}


/*
    Initializes everthing related to the flocks
        >also initializes the boids linked list for each flock
 */
void InitFlock(t_jit_boids3d *flockPtr)
{
    //General initialization
    flockPtr->number            = kNumBoids*MAX_FLOCKS;	//added init for jitter object
    flockPtr->neighbors			= kNumNeighbors;
    
    flockPtr->flyrect[top]		= kFlyRectTop;
    flockPtr->flyrect[left]		= kFlyRectLeft;
    flockPtr->flyrect[bottom]	= kFlyRectBottom;
    flockPtr->flyrect[right]	= kFlyRectRight;
    flockPtr->flyrect[front]	= kFlyRectFront;
    flockPtr->flyrect[back]		= kFlyRectBack;
    
    flockPtr->attractorLL = NULL;
    flockPtr->numAttractors = 0;
    flockPtr->allowNeighborsFromDiffFlock = 1;
    
    
    //Flock specific initialization
    for(int i=0; i<MAX_FLOCKS; i++){
        
        //create the linked list
        BoidPtr newLL = InitLL(flockPtr, kNumBoids, i);
        //TODO: error checking here breaks the external
        
        flockPtr->flockLL[i] = newLL; //add LL to the array
        
        //default values
        flockPtr->minspeed[i]			= kMinSpeed;
        flockPtr->maxspeed[i]			= kMaxSpeed;
        flockPtr->center[i]             = kCenterWeight;
        flockPtr->attract[i]			= kAttractWeight;
        flockPtr->match[i]				= kMatchWeight;
        flockPtr->avoid[i]				= kAvoidWeight;
        flockPtr->repel[i]				= kWallsWeight;
        flockPtr->edgedist[i]			= kEdgeDist;
        flockPtr->speed[i]				= kDefaultSpeed;
        flockPtr->inertia[i]			= kInertiaFactor;
        flockPtr->accel[i]              = kAccelFactor;
        flockPtr->neighborRadius[i]     = kNRadius;
    }
}


/*
 Initializes a LL. Returns a pointer to its head
 */
BoidPtr InitLL(t_jit_boids3d *flockPtr, long numBoids, int flockID)
{
    BoidPtr head = NULL;
    for(int i=0; i < numBoids; i++){
        
        BoidPtr theBoid = InitBoid(flockPtr); //make a new boid
        if(!theBoid){
            printf("ERROR: failed to malloc a boid\n");
            exit(1);
        }
        
        //Add the boid to the LL or make it the head (if nothing has already been added)
        if (head == NULL){
            head = theBoid;
            theBoid->nextBoid = NULL;
        }else{
            theBoid->nextBoid = head;
            head = theBoid;
        }
        
        //update number of boids in the flock
        flockPtr->boidCount[flockID]++;
        
    }
    
    return head;
}


/*
    Allocates memory for a Boid struct and returns a pointer to it
 */
BoidPtr InitBoid(t_jit_boids3d *flockPtr)
{
    struct Boid * theBoid = (struct Boid *)malloc(sizeof(struct Boid));
    
    if(!theBoid){
        return NULL;
    }
    
    theBoid->age = 0; //set age to 0
    
    //initialize struct variables
    theBoid->oldPos[x] = 0.0;
    theBoid->oldPos[y] = 0.0;
    theBoid->oldPos[z] = 0.0;
    
    theBoid->newPos[x] = 0.0;
    theBoid->newPos[y] = 0.0;
    theBoid->newPos[z] = 0.0;
    
    theBoid->oldDir[x] = 0.0;
    theBoid->oldDir[y] = 0.0;
    theBoid->oldDir[z] = 0.0;
    
    theBoid->newDir[x] = 0.0;
    theBoid->newDir[y] = 0.0;
    theBoid->newDir[z] = 0.0;
    
    theBoid->speed = 0.0;
    
    theBoid->newPos[x] = theBoid->oldPos[x] = (kFlyRectScalingFactor*RandomInt(flockPtr->flyrect[right],flockPtr->flyrect[left]));		// set random location within flyrect
    theBoid->newPos[y] = theBoid->oldPos[y] = (kFlyRectScalingFactor*RandomInt(flockPtr->flyrect[bottom], flockPtr->flyrect[top]));
    theBoid->newPos[z] = theBoid->oldPos[z] = (kFlyRectScalingFactor*RandomInt(flockPtr->flyrect[back], flockPtr->flyrect[front]));
    double rndAngle = RandomInt(0, 360) * flockPtr->d2r;		// set velocity from random angle
    theBoid->newDir[x] = jit_math_sin(rndAngle);
    theBoid->newDir[y] = jit_math_cos(rndAngle);
    theBoid->newDir[z] = (jit_math_cos(rndAngle) + jit_math_sin(rndAngle)) * 0.5;
    theBoid->speed = (kMaxSpeed + kMinSpeed) * 0.5;
    
    for(int j=0; j<kMaxNeighbors;j++) {
        theBoid->neighbor[j] = 0;
        theBoid->neighborDistSqr[j] = 0.0;
    }
    
    return theBoid;
}


/*
    Mallocs memory for an attractor and returns a pointer to it
 */
AttractorPtr InitAttractor(t_jit_boids3d *flockPtr)
{
    struct Attractor * theAttractor = (struct Attractor *)malloc(sizeof(struct Attractor));
    
    if(!theAttractor){
        return NULL;
    }
    
    //TODO: make the attractor at a specific point (not always the origin)
    theAttractor->loc[0] = 0.0;
    theAttractor->loc[1] = 0.0;
    theAttractor->loc[2] = 0.0;
    
    return theAttractor;
}


/*
    Calculates the total number of boids across all flocks
 */
int CalcNumBoids(t_jit_boids3d *flockPtr)
{
    int boidsCounter = 0;
    for (int i=0; i<MAX_FLOCKS; i++){
        boidsCounter += flockPtr->boidCount[i];
    }
    return boidsCounter;
}


/*
    Do initialization on startup of the jitter object
 */
t_jit_boids3d *jit_boids3d_new(void)
{
	t_jit_boids3d *flockPtr;
    
	if ((flockPtr=(t_jit_boids3d *)jit_object_alloc(_jit_boids3d_class))) {
		
		flockPtr->flyRectCount		= 6;
		//flockPtr->attractPtCount	= 1;
		flockPtr->mode	 			= 0;
		
		//init boids params
		InitFlock(flockPtr);
		
		flockPtr->d2r = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068/180.0;
		flockPtr->r2d = 180.0/3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068;
	} else {
		flockPtr = NULL;
	}
	return flockPtr;
}


/*
 Free the linked lists containing all the boids
    >NOTE: this replaces the previous method: jit_boids3d_free
 */
void freeFlocks(t_jit_boids3d *flockPtr)
{
    for(int i=0; i<MAX_FLOCKS; i++){ //we're clearing each flock
        
        if(flockPtr->flockLL[i] == NULL){ //ensure that this flock is populated
            continue;
        }
        
        BoidPtr iterator = flockPtr->flockLL[i];
        BoidPtr deletor = iterator;
        
        do{ //traverse the LL and free the memory for each boid
            iterator = iterator->nextBoid;
            
            free(deletor);
            deletor = iterator;
            
        }while (iterator);
        
        //mark the head null
        flockPtr->flockLL[i] = NULL;
        
    }
}
