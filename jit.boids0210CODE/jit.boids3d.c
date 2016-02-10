/*
 jit.boids3d.c  12/15/2005 wesley smith
 
 adapted from:
 boids3d 08/2005 a.sier / jasch adapted from boids by eric singer ÔøΩ 1995-2003 eric l. singer
 free for non-commercial use
 
 modified 070207 by sier to set boids pos, etc
 
 Modified by Jack Truskowski and Grace Handler 2015 to add functionality for multiple flocks
 -- commented methods that were heavily edited
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
const int           kBoidMaxAge     = 100;
const long			kNumBoids		= 0;	// number of boids
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
    
} Boid, *BoidPtr;

typedef struct _jit_boids3d
{
	t_object		ob;
	char			mode;
	long			number;
	long			neighbors;
	
	double			flyrect[6];
	long			flyRectCount;
    
    int             allowNeighborsFromDiffFlock;
	
    int             boidCount[MAX_FLOCKS];
    int             flockID[MAX_FLOCKS];
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
	
    //XYZ for each flock
    
    double          tempCenterPt[3];
	long			centerPtCount;
	
	double			attractpt[3];
	long			attractPtCount;
	
	BoidPtr			boid;
	
	double 			d2r;
	double			r2d;
    
} t_jit_boids3d;

void *_jit_boids3d_class;

t_jit_err jit_boids3d_init(void);
t_jit_boids3d *jit_boids3d_new(void);
void jit_boids3d_free(t_jit_boids3d *flockPtr);
t_jit_err jit_boids3d_matrix_calc(t_jit_boids3d *flockPtr, void *inputs, void *outputs);

void jit_boids3d_calculate_ndim(t_jit_boids3d *flockPtr, long dimcount, long *dim, long planecount,
                                t_jit_matrix_info *out_minfo, char *bop);

//attribute methods
t_jit_err jit_boids3d_neighbors(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
t_jit_err jit_boids3d_minspeed(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
t_jit_err jit_boids3d_nradius(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
t_jit_err jit_boids3d_inertia(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
t_jit_err jit_boids3d_number(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
//t_jit_err jit_boids3d_initialNumber(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
//t_jit_err jit_boids3d_boidCount(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
t_jit_err jit_boids3d_maxspeed(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
t_jit_err jit_boids3d_center(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
t_jit_err jit_boids3d_attract(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
t_jit_err jit_boids3d_match(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
t_jit_err jit_boids3d_avoid(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
t_jit_err jit_boids3d_repel(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
t_jit_err jit_boids3d_edgedist(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv); //NOT USED
t_jit_err jit_boids3d_speed(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);
t_jit_err jit_boids3d_accel(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv);

//Boids specific methods
void InitFlock(t_jit_boids3d *flockPtr);
void Flock_donumBoids(t_jit_boids3d *flockPtr, long numBoids);
void Flock_resetBoids(t_jit_boids3d *flockPtr);
void Flock_resetNewBoids(t_jit_boids3d *flockPtr, int startingIndex);

void FlightStep(t_jit_boids3d *flockPtr);
void FindFlockCenter(t_jit_boids3d *flockPtr, long theBoid);
float MatchAndAvoidNeighbors(t_jit_boids3d *flockPtr, long theBoid, double *matchNeighborVel, double *avoidNeighborVel);
void SeekPoint(t_jit_boids3d *flockPtr, long theBoid, double *seekPt, double* seekDir);
void AvoidWalls(t_jit_boids3d *flockPtr, long theBoid, double *wallVel);
char InFront(BoidPtr theBoid, BoidPtr neighbor);
void NormalizeVelocity(double *direction);
double RandomInt(double minRange, double maxRange);
double DistSqrToPt(double *firstPoint, double *secondPoint);


t_jit_err jit_boids3d_init(void)
{
    long attrflags=0;
	t_jit_object *attr,*mop,*o;
	t_symbol *atsym;
	
	atsym = gensym("jit_attr_offset");
	
	_jit_boids3d_class = jit_class_new("jit_boids3d",(method)jit_boids3d_new,(method)jit_boids3d_free,
                                       sizeof(t_jit_boids3d),0L);
    
	//add mop
	mop = jit_object_new(_jit_sym_jit_mop,0,1);
	o = jit_object_method(mop,_jit_sym_getoutput,1);
	jit_attr_setlong(o,_jit_sym_dimlink,0);
	
	jit_class_addadornment(_jit_boids3d_class,mop);
	//add methods
	jit_class_addmethod(_jit_boids3d_class, (method)jit_boids3d_matrix_calc, 		"matrix_calc", 		A_CANT, 0L);
	jit_class_addmethod(_jit_boids3d_class, (method)Flock_resetBoids, 				"reset", 			A_USURP_LOW, 0L);
    
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
    
	//attractpt.
	attr = jit_object_new(_jit_sym_jit_attr_offset_array,"attractpt",_jit_sym_float64,4,attrflags,
                          (method)0L,(method)0L,calcoffset(t_jit_boids3d,attractPtCount),calcoffset(t_jit_boids3d,attractpt));
	jit_class_addattr(_jit_boids3d_class,attr);
	
	jit_class_register(_jit_boids3d_class);
    
	return JIT_ERR_NONE;
    
}

//boids attribute methods

t_jit_err jit_boids3d_neighbors(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
	flockPtr->neighbors = (double)MIN(jit_atom_getfloat(argv), kMaxNeighbors);
}

t_jit_err jit_boids3d_nradius(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
    int flockID = (int)jit_atom_getfloat(argv+1);
	flockPtr->neighborRadius[flockID] = (double)MAX(jit_atom_getfloat(argv), 0.0);
}

t_jit_err jit_boids3d_minspeed(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
    int flockID = (int)jit_atom_getfloat(argv+1);
	flockPtr->minspeed[flockID] = (double)MAX(jit_atom_getfloat(argv), 0.000001);
}

t_jit_err jit_boids3d_maxspeed(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
    int flockID = (int)jit_atom_getfloat(argv+1);
	flockPtr->maxspeed[flockID] = (double)MAX(jit_atom_getfloat(argv), 0.000001);
}

t_jit_err jit_boids3d_center(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
    int flockID = (int)jit_atom_getfloat(argv+1);
	flockPtr->center[flockID] = (double)MAX(jit_atom_getfloat(argv), 0.000001);
}

t_jit_err jit_boids3d_attract(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
    int flockID = (int)jit_atom_getfloat(argv+1);
	flockPtr->attract[flockID] = (double)MAX(jit_atom_getfloat(argv), 0.000001);
}

t_jit_err jit_boids3d_match(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
    int flockID = (int)jit_atom_getfloat(argv+1);
	flockPtr->match[flockID] = (double)MAX(jit_atom_getfloat(argv), 0.000001);
}

t_jit_err jit_boids3d_avoid(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
    int flockID = (int)jit_atom_getfloat(argv+1);
	flockPtr->avoid[flockID] = (double)MAX(jit_atom_getfloat(argv), 0.000001);
}

t_jit_err jit_boids3d_repel(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
    int flockID = (int)jit_atom_getfloat(argv+1);
	flockPtr->repel[flockID] = (double)MAX(jit_atom_getfloat(argv), 0.000001);
}

t_jit_err jit_boids3d_edgedist(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
    int flockID = (int)jit_atom_getfloat(argv+1);
	flockPtr->edgedist[flockID] = (double)MAX(jit_atom_getfloat(argv), 0.000001);
}

t_jit_err jit_boids3d_speed(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
    int flockID = (int)jit_atom_getfloat(argv+1);
	flockPtr->speed[flockID] = (double)MAX(jit_atom_getfloat(argv), 0.000001);
}

t_jit_err jit_boids3d_inertia(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
	double val = (double)jit_atom_getfloat(argv);
    int flockID =(int)jit_atom_getfloat(argv+1);
	
	if(val == 0.0)
		flockPtr->inertia[flockID] = 0.000001;
	else
		flockPtr->inertia[flockID] = val;
}

t_jit_err jit_boids3d_accel(t_jit_boids3d *flockPtr, void *attr, long argc, t_atom *argv)
{
    int flockID = (int)jit_atom_getfloat(argv+1);
	flockPtr->accel[flockID] = (double)MAX(jit_atom_getfloat(argv), 0.000001);
}

/*
    This function is called when the number of boids in the maxpatch is changed. Creates and initializes a new array containing all boids
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
    
    //create new array for all boids in all flocks
    Boid *boidArray = (Boid *) jit_getbytes(sizeof(Boid)*newNumBoids);
    
    int index = 0;
    
    //for every boid, check if it needs to be deleted
    for (int i=0; i<flockPtr->number; i++){
        if (boidChanges[flockPtr->boid[i].flockID] < 0){ //not adding this boid, increment boidChanges by 1
            //we're deleting this boid
            boidChanges[flockPtr->boid[i].flockID]++;
        }
        else{
            //add this boid to the new array
            boidArray[index] = flockPtr->boid[i]; //add this boid to the new total boids array
            index++;
        }
    }
    
    //index of the last boid in boidArray
    int indexOfNewBoids = index;
    
    //assigns flock IDs to new boids
    for (int i=0; i<6; i++){
        while (boidChanges[i] > 0){
            //add a new boid
            boidArray[index].flockID = i;
            boidChanges[i]--;
            index++;
        }
    }
    
    //THIS IS WHERE WE CAN QUICKSORT THE BOID ARRAY IN THE CODE SINCE BUBBLE SORT (IN MAX) IS INEFFICIENT
    
    //Assign the new array to be flockPtr->boid
    flockPtr->number = newNumBoids;
    Flock_donumBoids(flockPtr, newNumBoids);
    flockPtr->boid = boidArray;
    
    
    Flock_resetNewBoids(flockPtr, indexOfNewBoids);
}


/*
    Prepares the output matrix and sends it back to the max patch
 */
t_jit_err jit_boids3d_matrix_calc(t_jit_boids3d *flockPtr, void *inputs, void *outputs)
{
	t_jit_err err=JIT_ERR_NONE;
	long out_savelock;
	t_jit_matrix_info out_minfo;
	char *out_bp;
	long i,dimcount,planecount,dim[JIT_MATRIX_MAX_DIMCOUNT];
	void *out_matrix;
	
	out_matrix 	= jit_object_method(outputs,_jit_sym_getindex,0);
    
	if (flockPtr&&out_matrix) {
		out_savelock = (long) jit_object_method(out_matrix,_jit_sym_lock,1);
		
		jit_object_method(out_matrix,_jit_sym_getinfo,&out_minfo);
		
        //dimensions of the output matrix (number of boids x 1)
		out_minfo.dim[0] = flockPtr->number;
		out_minfo.dim[1] = 1;
		out_minfo.type = _jit_sym_float32; //outputting floating point numbers
		
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
		
		jit_object_method(out_matrix,_jit_sym_setinfo,&out_minfo);
		jit_object_method(out_matrix,_jit_sym_getinfo,&out_minfo);
		
		jit_object_method(out_matrix,_jit_sym_getdata,&out_bp);
		
		if (!out_bp) { err=JIT_ERR_INVALID_OUTPUT; goto out;}
		
		//get dimensions/planecount
		dimcount   = out_minfo.dimcount;
		planecount = out_minfo.planecount;
		
		for (i=0;i<dimcount;i++) {
			dim[i] = out_minfo.dim[i];
		}
        
		jit_boids3d_calculate_ndim(flockPtr, dimcount, dim, planecount, &out_minfo, out_bp);
	} else {
		return JIT_ERR_INVALID_PTR;
	}
	
out: //output the matrix
	jit_object_method(out_matrix,gensym("lock"),out_savelock);
	return err;
}

//GRACE AND JACK CHANGED THIS HEAVILY
void jit_boids3d_calculate_ndim(t_jit_boids3d *flockPtr, long dimcount, long *dim, long planecount,
                                t_jit_matrix_info *out_minfo, char *bop)
{
	long i, k;
	float *fop;
	BoidPtr boid;
	double 	tempNew_x, tempNew_y, tempNew_z;
	double 	tempOld_x, tempOld_y, tempOld_z;
	double	delta_x, delta_y, delta_z, azi, ele, speed;
	
    //do a step in the simulation
	FlightStep(flockPtr);
	
	//copy pointer so we don't have to dereference in the for loop
	boid = flockPtr->boid;
	fop = (float *)bop; //contains the planes
	
    //apply the correct mode and fill the planes with the correct info
    switch(flockPtr->mode) { // newpos
		case 0:
			for(i=0; i < dim[0]; i++) {
				fop[0] = boid[i].newPos[x];
				fop[1] = boid[i].newPos[y];
				fop[2] = boid[i].newPos[z];
                fop[3] = boid[i].flockID;
				
				fop += planecount;
			}
            break;
		case 1: //newpos + oldpos
			for(i=0; i < dim[0]; i++) {
				fop[0] = boid[i].newPos[x];
				fop[1] = boid[i].newPos[y];
				fop[2] = boid[i].newPos[z];
                fop[3] = boid[i].flockID;
				fop[4] = boid[i].oldPos[x];
				fop[5] = boid[i].oldPos[y];
				fop[6] = boid[i].oldPos[z];
				
				fop += planecount;
			}
            break;
		case 2: //newpos + oldpos + speed-azimuth-elevation
			for(i=0; i < dim[0]; i++) {
				tempNew_x = boid[i].newPos[x];
				tempNew_y = boid[i].newPos[y];
				tempNew_z = boid[i].newPos[z];
				tempOld_x = boid[i].oldPos[x];
				tempOld_y = boid[i].oldPos[y];
				tempOld_z = boid[i].oldPos[z];
				
				delta_x = tempNew_x - tempOld_x;
				delta_y = tempNew_y - tempOld_y;
				delta_z = tempNew_z - tempOld_z;
				azi = jit_math_atan2(delta_z, delta_x) * flockPtr->r2d;
				ele = jit_math_atan2(delta_y, delta_x) * flockPtr->r2d;
				speed = jit_math_sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);
				
				fop[0] = tempNew_x;
				fop[1] = tempNew_y;
				fop[2] = tempNew_z;
                fop[3] = boid[i].flockID;
				fop[4] = tempOld_x;
				fop[5] = tempOld_y;
				fop[6] = tempOld_z;
				fop[7] = speed;
				fop[8] = azi;
				fop[9] = ele;
				
				fop += planecount;
			}
            break;
	}

}

/*
    Calculates each boid's new velocity and updates position
 */
void FlightStep(t_jit_boids3d *flockPtr)
{
	double			goCenterVel[3];
	double			goAttractVel[3];
	double			matchNeighborVel[3];
    double			avoidWallsVel[3];
	double			avoidNeighborVel[3];
	double			avoidNeighborSpeed;
	const double	zeroVel[3]	= {0.0, 0.0, 0.0};
	long			i;
    
	for (i = 0; i <  flockPtr->number; i++) {						// save position and velocity
		flockPtr->boid[i].oldPos[x] = flockPtr->boid[i].newPos[x];
		flockPtr->boid[i].oldPos[y] = flockPtr->boid[i].newPos[y];
		flockPtr->boid[i].oldPos[z] = flockPtr->boid[i].newPos[z];
		 
		flockPtr->boid[i].oldDir[x] = flockPtr->boid[i].newDir[x];
		flockPtr->boid[i].oldDir[y] = flockPtr->boid[i].newDir[y];
		flockPtr->boid[i].oldDir[z] = flockPtr->boid[i].newDir[z];
	}
    
    //for every boid in the simulation...
	for (i = 0; i < flockPtr->number; i++) {
        
        int flockID = flockPtr->boid[i].flockID;
        
        FindFlockCenter(flockPtr, i); //for centering instinct
        
		if (flockPtr->neighbors > 0) {							// get all velocity components
			avoidNeighborSpeed = MatchAndAvoidNeighbors(flockPtr, i, matchNeighborVel,  avoidNeighborVel);
		} else {
			matchNeighborVel[x] = zeroVel[x];
			matchNeighborVel[y] = zeroVel[y];
			matchNeighborVel[z] = zeroVel[z];
            
			avoidNeighborVel[x] = zeroVel[x];
			avoidNeighborVel[y] = zeroVel[y];
			avoidNeighborVel[z] = zeroVel[z];
			
			avoidNeighborSpeed = 0;
		}
        
        //update velocity to include centering and attracting instincts
		SeekPoint(flockPtr, i, flockPtr->tempCenterPt, goCenterVel);
		SeekPoint(flockPtr, i, flockPtr->attractpt, goAttractVel);
        
        //AVOID WALLS USED TO BE DONE HERE AS A COMPONENT OF VELOCITY
		//AvoidWalls(flockPtr, i, avoidWallsVel);
        
        
		// compute resultant velocity using weights and inertia
		flockPtr->boid[i].newDir[x] = flockPtr->inertia[flockID] * (flockPtr->boid[i].oldDir[x]) +
        (flockPtr->center[flockID] * goCenterVel[x] +
         flockPtr->attract[flockID] * goAttractVel[x] +
         flockPtr->match[flockID] * matchNeighborVel[x] +
         flockPtr->avoid[flockID] * avoidNeighborVel[x]) / flockPtr->inertia[flockID];
		flockPtr->boid[i].newDir[y] = flockPtr->inertia[flockID] * (flockPtr->boid[i].oldDir[y]) +
        (flockPtr->center[flockID] * goCenterVel[y] +
         flockPtr->attract[flockID] * goAttractVel[y] +
         flockPtr->match[flockID] * matchNeighborVel[y] +
         flockPtr->avoid[flockID] * avoidNeighborVel[y]) / flockPtr->inertia[flockID];
		flockPtr->boid[i].newDir[z] = flockPtr->inertia[flockID] * (flockPtr->boid[i].oldDir[z]) +
        (flockPtr->center[flockID] * goCenterVel[z] +
         flockPtr->attract[flockID] * goAttractVel[z] +
         flockPtr->match[flockID] * matchNeighborVel[z] +
         flockPtr->avoid[flockID] * avoidNeighborVel[z]) / flockPtr->inertia[flockID];
		NormalizeVelocity(flockPtr->boid[i].newDir);	// normalize velocity so its length is unity
        
		// set to avoidNeighborSpeed bounded by minspeed and maxspeed
		if ((avoidNeighborSpeed >= flockPtr->minspeed[flockID]) &&
            (avoidNeighborSpeed <= flockPtr->maxspeed[flockID]))
			flockPtr->boid[i].speed = avoidNeighborSpeed;
		else if (avoidNeighborSpeed > flockPtr->maxspeed[flockID])
			flockPtr->boid[i].speed = flockPtr->maxspeed[flockID];
		else
			flockPtr->boid[i].speed = flockPtr->minspeed[flockID];
        
        //bounce back from walls if the boid is beyond the limit of the flyrect
        AvoidWalls(flockPtr, i, flockPtr->boid[i].newDir);
        
		// calculate new position, applying speed
		flockPtr->boid[i].newPos[x] += flockPtr->boid[i].newDir[x] * flockPtr->boid[i].speed * (flockPtr->speed[flockID] / 100.0);
		flockPtr->boid[i].newPos[y] += flockPtr->boid[i].newDir[y] * flockPtr->boid[i].speed * (flockPtr->speed[flockID] / 100.0);
		flockPtr->boid[i].newPos[z] += flockPtr->boid[i].newDir[z] * flockPtr->boid[i].speed * (flockPtr->speed[flockID] / 100.0);
        
	}
}

/*
    Calculates the center of a flock, saves it in flockPtr->centerPt
    TODO: is this being done for every boid in a flock? Should be done only once
 */
void FindFlockCenter(t_jit_boids3d *flockPtr, long theBoid)
{
    int flockID = flockPtr->boid[theBoid].flockID;
    double totalH = 0, totalV = 0, totalD = 0;
    int neighborsCount = 0;
    
    for(int i=0; i<flockPtr->number; i++){
        
        //grab the boid and calculate distance
        Boid currBoid = flockPtr->boid[i];
        double dist = sqrt(DistSqrToPt(flockPtr->boid[theBoid].oldPos, flockPtr->boid[i].oldPos));
        
        //ensure that the boid is not looking at itself and we're not overfilling neighborList
        if(dist < flockPtr->neighborRadius[flockID] && dist > 0.0 && neighborsCount < kMaxNeighbors){
            
            //check to ensure this boid is allowed / in same flock
            if (flockPtr->allowNeighborsFromDiffFlock == 0 && currBoid.flockID != flockID){
                continue;
            }
            
            //this boid is a neighbor
            neighborsCount++;
            totalH += flockPtr->boid[i].oldPos[x];
            totalV += flockPtr->boid[i].oldPos[y];
            totalD += flockPtr->boid[i].oldPos[z];
        }
        
    }
    
    if(neighborsCount > 0){ //get the average position of all boids in the flock
        flockPtr->tempCenterPt[x] = (double)	(totalH / neighborsCount);
        flockPtr->tempCenterPt[y] = (double)	(totalV / neighborsCount);
        flockPtr->tempCenterPt[z] = (double)	(totalD / neighborsCount);
    }else{ //only boid in flock, its position is the center point
        flockPtr->tempCenterPt[x] = flockPtr->boid[theBoid].oldPos[x];
        flockPtr->tempCenterPt[y] = flockPtr->boid[theBoid].oldPos[y];
        flockPtr->tempCenterPt[z] = flockPtr->boid[theBoid].oldPos[z];
    }

}

/*
    Computes the matching and avoidance of neighbors instinct
 */
float MatchAndAvoidNeighbors(t_jit_boids3d *flockPtr, long theBoid, double *matchNeighborVel, double *avoidNeighborVel)
{
	long			i, j;
    Boid            neighbor;
	double			distSqr;
	double			dist, distH, distV,distD;
	double			tempSpeed;
	short			numClose = 0;
	double			totalVel[3] = {0.0,0.0,0.0};
    int flockID = flockPtr->boid[theBoid].flockID;
    
    Boid neighborList[kMaxNeighbors];
    int neighborsCount = 0; //counter to keep track of how many neighbors we've found
    
	/**********************/
	/* Find the neighbors */
	/**********************/
    
    //loop through all boids
    for(i=0; i<flockPtr->number; i++){
        
        //grab the boid and calculate distance
        Boid currBoid = flockPtr->boid[i];
        dist = sqrt(DistSqrToPt(flockPtr->boid[theBoid].oldPos, flockPtr->boid[i].oldPos));
        
        //ensure that the boid is not looking at itself and we're not overfilling neighborList
        if(dist < flockPtr->neighborRadius[flockID] && dist > 0.0 && neighborsCount < kMaxNeighbors){
            
            //check to ensure this boid is allowed / in same flock
            if (flockPtr->allowNeighborsFromDiffFlock == 0 && currBoid.flockID != flockID){
                continue;
            }
            
            //this boid is a neighbor
            neighborList[neighborsCount] = currBoid;
            neighborsCount++;
        }
        
    }
    
    
    
    
//	/* special case of one neighbor */
//	if (flockPtr->neighbors == 1) {
//		flockPtr->boid[theBoid].neighborDistSqr[0] = kMaxLong;
//        
//		for (i = 0; i < flockPtr->number; i++) {
//			if (i != theBoid) {
//				distSqr = DistSqrToPt(flockPtr->boid[theBoid].oldPos, flockPtr->boid[i].oldPos);
//				
//				/* if this one is closer than the closest so far, then remember it */
//				if (flockPtr->boid[theBoid].neighborDistSqr[0] > distSqr) {
//					flockPtr->boid[theBoid].neighborDistSqr[0] = distSqr;
//					flockPtr->boid[theBoid].neighbor[0] = i;
//				}
//			}
//		}
//	}
//	/* more than one neighbor */
//	else {
//		for (j = 0; j < flockPtr->neighbors; j++)
//			flockPtr->boid[theBoid].neighborDistSqr[j] = kMaxLong;
//		
//		for (i = 0 ; i < flockPtr->number; i++) {
//			/* if this one is not me... */
//			if (i != theBoid) {
//				distSqr = DistSqrToPt(flockPtr->boid[theBoid].oldPos, flockPtr->boid[i].oldPos);
//                
//				/* if distSqr is less than the distance at the bottom of the array, sort into array */
//				if (distSqr < flockPtr->boid[theBoid].neighborDistSqr[flockPtr->neighbors-1]) {
//					j = flockPtr->neighbors - 1;
//                    
//					/* sort distSqr in to keep array in size order, smallest first */
//					while ((distSqr < flockPtr->boid[theBoid].neighborDistSqr[j-1]) && (j > 0)) {
//						flockPtr->boid[theBoid].neighborDistSqr[j] = flockPtr->boid[theBoid].neighborDistSqr[j - 1];
//						flockPtr->boid[theBoid].neighbor[j] = flockPtr->boid[theBoid].neighbor[j - 1];
//						j--;
//					}
//					flockPtr->boid[theBoid].neighborDistSqr[j] = distSqr;
//					flockPtr->boid[theBoid].neighbor[j] = i;
//				}
//			}
//		}
//	}
    
	/*********************************/
	/* Match and avoid the neighbors */
	/*********************************/
    
	matchNeighborVel[x] = 0;
	matchNeighborVel[y] = 0;
	matchNeighborVel[z] = 0;
	
	// set tempSpeed to old speed
	tempSpeed = flockPtr->boid[theBoid].speed;
	
	for (i = 0; i < neighborsCount; i++) {
		neighbor = neighborList[i];
		
		// calculate matchNeighborVel by averaging the neighbor velocities
		matchNeighborVel[x] += neighbor.oldDir[x];
		matchNeighborVel[y] += neighbor.oldDir[y];
		matchNeighborVel[z] += neighbor.oldDir[z];
        
        distH = neighbor.oldPos[x] - flockPtr->boid[theBoid].oldPos[x];
        distV = neighbor.oldPos[y] - flockPtr->boid[theBoid].oldPos[y];
        distD = neighbor.oldPos[z] - flockPtr->boid[theBoid].oldPos[z];
        
        if(dist == 0.0) dist = 0.0000001;
        
        /*  This can be used to give the boids a preferred distance from their neighbors
        totalVel[x] = totalVel[x] - distH - (distH * ((float) flockPtr->neighborRadius[flockID] / (dist)));
        totalVel[y] = totalVel[y] - distV - (distV * ((float) flockPtr->neighborRadius[flockID] / (dist)));
        totalVel[z] = totalVel[z] - distD - (distV * ((float) flockPtr->neighborRadius[flockID] / (dist)));
         */
        
        numClose++;
		
		
		if (InFront(&(flockPtr->boid[theBoid]), &(neighbor))) {	// adjust speed
				tempSpeed /= (flockPtr->accel[flockID] / 100.0);
		}
		else {
				tempSpeed *= (flockPtr->accel[flockID] / 100.0);
		}
	}
	if (numClose) {
		avoidNeighborVel[x] = totalVel[x] / numClose;
		avoidNeighborVel[y] = totalVel[y] / numClose;
		avoidNeighborVel[z] = totalVel[z] / numClose;
		NormalizeVelocity(matchNeighborVel);
	}
	else {
		avoidNeighborVel[x] = 0;
		avoidNeighborVel[y] = 0;
		avoidNeighborVel[z] = 0;
	}
	return(tempSpeed);
}

//GRACE AND JACK CHANGED THIS HEAVILY
//modified to take the returned velocity as an argument instead
void SeekPoint(t_jit_boids3d *flockPtr, long theBoid, double *seekPt, double* seekDir)
{
	seekDir[x] = seekPt[x] - flockPtr->boid[theBoid].oldPos[x];
	seekDir[y] = seekPt[y] - flockPtr->boid[theBoid].oldPos[y];
	seekDir[z] = seekPt[z] - flockPtr->boid[theBoid].oldPos[z];
	NormalizeVelocity(seekDir);
}

void AvoidWalls(t_jit_boids3d *flockPtr, long theBoid, double *wallVel)
{
	double		testPoint[3];
    
//	wallVel[x] = 0.0;
//	wallVel[y] = 0.0;
//	wallVel[z] = 0.0;
    
	/* calculate test point in front of the nose of the boid */
	/* distance depends on the boid's speed and the avoid edge constant */
	testPoint[x] = flockPtr->boid[theBoid].oldPos[x] + flockPtr->boid[theBoid].newDir[x] * (flockPtr->boid[theBoid].speed * (flockPtr->speed[flockPtr->boid[theBoid].flockID] / 100.0));// * flockPtr->edgedist[flockPtr->boid[theBoid].flockID];
    testPoint[y] = flockPtr->boid[theBoid].oldPos[y] + flockPtr->boid[theBoid].newDir[y] * (flockPtr->boid[theBoid].speed * (flockPtr->speed[flockPtr->boid[theBoid].flockID] / 100.0));// * flockPtr->edgedist[flockPtr->boid[theBoid].flockID];fl
    testPoint[z] = flockPtr->boid[theBoid].oldPos[z] + flockPtr->boid[theBoid].newDir[z] * (flockPtr->boid[theBoid].speed * (flockPtr->speed[flockPtr->boid[theBoid].flockID] / 100.0));// * flockPtr->edgedist[flockPtr->boid[theBoid].flockID];

    
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

//GRACE AND JACK CHANGED THIS HEAVILY
void Flock_resetNewBoids(t_jit_boids3d *flockPtr, int startingIndex)
{
    long i, j;
	double rndAngle;
	
	for (i = startingIndex; i <  flockPtr->number; i++) { // init everything to 0.0
		flockPtr->boid[i].oldPos[x] = 0.0;
		flockPtr->boid[i].oldPos[y] = 0.0;
		flockPtr->boid[i].oldPos[z] = 0.0;
        
		flockPtr->boid[i].newPos[x] = 0.0;
		flockPtr->boid[i].newPos[y] = 0.0;
		flockPtr->boid[i].newPos[z] = 0.0;
		
		flockPtr->boid[i].oldDir[x] = 0.0;
		flockPtr->boid[i].oldDir[y] = 0.0;
		flockPtr->boid[i].oldDir[z] = 0.0;
		
		flockPtr->boid[i].newDir[x] = 0.0;
		flockPtr->boid[i].newDir[y] = 0.0;
		flockPtr->boid[i].newDir[z] = 0.0;
		
		flockPtr->boid[i].speed = 0.0;
		
		for(j=0; j<kMaxNeighbors;j++) {
			flockPtr->boid[i].neighbor[j] = 0;
			flockPtr->boid[i].neighborDistSqr[j] = 0.0;
		}
	}
	
	for (i = startingIndex; i <  flockPtr->number; i++) {				// set the initial locations and velocities of the boids
		flockPtr->boid[i].newPos[x] = flockPtr->boid[i].oldPos[x] = RandomInt(flockPtr->flyrect[right],flockPtr->flyrect[left]);		// set random location within flyrect
		flockPtr->boid[i].newPos[y] = flockPtr->boid[i].oldPos[y] = RandomInt(flockPtr->flyrect[bottom], flockPtr->flyrect[top]);
		flockPtr->boid[i].newPos[z] = flockPtr->boid[i].oldPos[z] = RandomInt(flockPtr->flyrect[back], flockPtr->flyrect[front]);
		rndAngle = RandomInt(0, 360) * flockPtr->d2r;		// set velocity from random angle
		flockPtr->boid[i].newDir[x] = jit_math_sin(rndAngle);
		flockPtr->boid[i].newDir[y] = jit_math_cos(rndAngle);
		flockPtr->boid[i].newDir[z] = (jit_math_cos(rndAngle) + jit_math_sin(rndAngle)) * 0.5;
		flockPtr->boid[i].speed = (kMaxSpeed + kMinSpeed) * 0.5;
	}
}

//UPDATED
void Flock_resetBoids(t_jit_boids3d *flockPtr)
{
	long i, j;
	double rndAngle;
	
	for (i = 0; i <  flockPtr->number; i++) { // init everything to 0.0
		flockPtr->boid[i].oldPos[x] = 0.0;
		flockPtr->boid[i].oldPos[y] = 0.0;
		flockPtr->boid[i].oldPos[z] = 0.0;
        
		flockPtr->boid[i].newPos[x] = 0.0;
		flockPtr->boid[i].newPos[y] = 0.0;
		flockPtr->boid[i].newPos[z] = 0.0;
		
		flockPtr->boid[i].oldDir[x] = 0.0;
		flockPtr->boid[i].oldDir[y] = 0.0;
		flockPtr->boid[i].oldDir[z] = 0.0;
		
		flockPtr->boid[i].newDir[x] = 0.0;
		flockPtr->boid[i].newDir[y] = 0.0;
		flockPtr->boid[i].newDir[z] = 0.0;
		
		flockPtr->boid[i].speed = 0.0;
		
		for(j=0; j<kMaxNeighbors;j++) {
			flockPtr->boid[i].neighbor[j] = 0;
			flockPtr->boid[i].neighborDistSqr[j] = 0.0;
		}
	}
	
	for (i = 0; i <  flockPtr->number; i++) {				// set the initial locations and velocities of the boids
		flockPtr->boid[i].newPos[x] = flockPtr->boid[i].oldPos[x] = RandomInt(flockPtr->flyrect[right],flockPtr->flyrect[left]);		// set random location within flyrect
		flockPtr->boid[i].newPos[y] = flockPtr->boid[i].oldPos[y] = RandomInt(flockPtr->flyrect[bottom], flockPtr->flyrect[top]);
		flockPtr->boid[i].newPos[z] = flockPtr->boid[i].oldPos[z] = RandomInt(flockPtr->flyrect[back], flockPtr->flyrect[front]);
		rndAngle = RandomInt(0, 360) * flockPtr->d2r;		// set velocity from random angle
		flockPtr->boid[i].newDir[x] = jit_math_sin(rndAngle);
		flockPtr->boid[i].newDir[y] = jit_math_cos(rndAngle);
		flockPtr->boid[i].newDir[z] = (jit_math_cos(rndAngle) + jit_math_sin(rndAngle)) * 0.5;
		flockPtr->boid[i].speed = (kMaxSpeed + kMinSpeed) * 0.5;
	}
}


//GRACE AND JACK CHANGED THIS HEAVILY
void InitFlock(t_jit_boids3d *flockPtr)
{
    //General initialization
    flockPtr->number            = kNumBoids;	//added init for jitter object
    flockPtr->neighbors			= kNumNeighbors;
    
    flockPtr->flyrect[top]		= kFlyRectTop;
    flockPtr->flyrect[left]		= kFlyRectLeft;
    flockPtr->flyrect[bottom]	= kFlyRectBottom;
    flockPtr->flyrect[right]	= kFlyRectRight;
    flockPtr->flyrect[front]	= kFlyRectFront;
    flockPtr->flyrect[back]		= kFlyRectBack;
    flockPtr->attractpt[x]		= (kFlyRectLeft + kFlyRectRight) * 0.5;
    flockPtr->attractpt[y]		= (kFlyRectTop + kFlyRectBottom) * 0.5;
    flockPtr->attractpt[z]		= (kFlyRectFront + kFlyRectBack) * 0.5;
    flockPtr->allowNeighborsFromDiffFlock = 1;
    
    Flock_donumBoids(flockPtr, flockPtr->number);
    
    //Set every boids flockID and age to 0 for initialization
    for (int i = 0; i <  flockPtr->number; i++){
        flockPtr->boid[i].flockID = 0;
        flockPtr->boid[i].age = 0;
    }

    
    //Flock specific initialization
    for(int i=0; i<MAX_FLOCKS; i++){
        flockPtr->minspeed[i]			= kMinSpeed;
        flockPtr->maxspeed[i]			= kMaxSpeed;
        flockPtr->center[i]             = kCenterWeight;
        flockPtr->attract[i]			= kAttractWeight;
        flockPtr->match[i]				= kMatchWeight;
        flockPtr->avoid[i]				= kAvoidWeight;
        flockPtr->repel[i]				= kWallsWeight;
        flockPtr->edgedist[i]			= kEdgeDist;
        flockPtr->speed[i]				= kSpeedupFactor;
        flockPtr->inertia[i]			= kInertiaFactor;
        flockPtr->accel[i]              = kAccelFactor;
        flockPtr->neighborRadius[i]     = kNRadius;
    }
    Flock_resetBoids(flockPtr);
}

//GRACE AND JACK CHANGED THIS HEAVILY
void Flock_donumBoids(t_jit_boids3d *flockPtr, long numBoids)
{
	//if not NULL, deallocate memory before assigning new memory
	//could check for whether flockPtr->number changed in value but not really
	//necessary as speed isn't the biggest issue in setting this param
	if(flockPtr->boid) {
		jit_boids3d_free(flockPtr);
	}
	
	flockPtr->boid = (Boid *) jit_getbytes(sizeof(Boid)*numBoids);
}


t_jit_boids3d *jit_boids3d_new(void)
{
	t_jit_boids3d *flockPtr;
    
	if ((flockPtr=(t_jit_boids3d *)jit_object_alloc(_jit_boids3d_class))) {
		
		flockPtr->flyRectCount		= 6;
		flockPtr->attractPtCount	= 3;
		flockPtr->mode	 			= 0;
		flockPtr->boid 				= NULL;
		
		//init boids params
		InitFlock(flockPtr);
		
		flockPtr->d2r = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068/180.0;
		flockPtr->r2d = 180.0/3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068;
	} else {
		flockPtr = NULL;
	}
	return flockPtr;
}

void jit_boids3d_free(t_jit_boids3d *flockPtr)
{
	//free bytes allocated for boids struct
	jit_freebytes((void *)flockPtr->boid, sizeof(Boid)*flockPtr->number);
}