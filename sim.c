/* Simulator for DTN under Random Waypoint mobility.
 * 
 * Nodes moves on a torus.
 * To change the FR, modify the function is_in_fr
 * To change the RR, modify the function is_better
 * The simulator will run and it will produce one log files:
 *  - log_stages.txt: log where each line represents a stages and its
 *    quantities needed to calculate the final metrics.
 *  - metrics.txt: different metrics as in the paper.
 *  - psi.txt: histogram of the limiting distribution function Psi.
 * 
 * If, for any reason you want to quit the simulation, hit Ctrl+c. The
 * log will be generated with the intermediate results. This is useful
 * to quit a simulation that takes too long due to the "infinite routing"
 * problem.
 *
 * Build with: 
 * 		gcc sim.c -Wall -lgsl -lgslcblas -lm -Ofast -o sim
 * 
 * author: riccardo.cavallari@unibo.it
 *
 * For the discretization of thetas:
 * use this formula to calculate the center of the bins
 * thetas=-pi+(pi/N)*(2*(1:1:N)-1);
 *
 * Changelog
 * 
 * 11.03.2018: added a parameter to use a constant seed for rand.
 *
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <unistd.h>
#include <sys/stat.h>
#include <gsl/gsl_histogram.h>

// Macro
#define TORUS_X			    100
#define TORUS_Y			    TORUS_X
#define ALLOC_MIN_SIZE  	10
#define NUM_ALLOC_STAGES  	100
#define HIST_BINS		    128
#define DIR_TAU             0.1     // to avoid infinite routing loop
#define min(a, b) (((a) < (b)) ? (a) : (b))

// Typedefs
struct Node {
	int nodeId;	    // id of the node
	double x;		// current x-position
	double y;		// current y-position
	double t;		// time to travel before next turn
	double d;		// direction of travel
	double dx;	    // x increment during a leg
	double dy;      // y increment during a leg
};

struct Packet {
	double d;		// current direction of travel
	double x;		// current x-position
	double y;		// current y-position
	int carrierId; 	// id of the current carrier node
};

typedef struct {
	double dir;		    // direction
	double tStart;	    // start time
	double tEnd;	    // end time
	double txDist;	    // distance covered by tx (can be zero)
	double txCost;      // tx cost (can be zero)
	unsigned int numHops;	// number of hops
}stage_t;

// Global varialbles
unsigned long int numOfNodes;
double lambda;  // density
double a;       // boundary function parameter
double b;       // boundary function parameter
double v0;      // nodes' speed
double r0;      // turning rate
double tw;      // parameter of f_D(x)
double simDur;  // duration of the simulation
double dT;      // duration of the simulation step
double dS;      // distance covered in dT
int    log_stages = 0; // 1 for logging stages, 0 for not
char   boundary;    // boundary function
char   potential;   // potential function
int stageOrTime;    // if 1 simDur is in number of stages; if 0, simDur is in time
int seedRand; // if -1 seed it time(NULL); if > -1, seed is seedRand 
FILE * f_stages;
FILE * f_metrics;
FILE * f_psi;
stage_t *stages = NULL;
unsigned int numAllocStages;	// sizeof(unsigned int) = 4
unsigned int numUsedStages;
gsl_histogram * hist = NULL;
	

/*
 * Generate samples of an exponentially distribute random variable
 * with a given rate.
 */
double next_time(double rate)
{	
    double r = 0;
	while (r == 0) {	
		r = (double)rand()/RAND_MAX;    // r is unif distr. in [0,1)
	}
	return -log(r)/rate;
}

/* 
 * Calculate the distance between two nodes on a torus.
 * sqrt(min(|x1 - x2|, w - |x1 - x2|)^2 + min(|y1 - y2|, h - |y1-y2|)^2)
 */
double dist_square(struct Node * n1, struct Node * n2)
{
	double absx = fabs(n1->x - n2->x);
	double absy = fabs(n1->y - n2->y);
	return (pow(min(absx, TORUS_X-absx),2) + pow(min(absy, TORUS_Y-absy),2));
}

/* 
 * Move nodes within the torus.
 */
void move_nodes(struct Node * nodeList)
{
  unsigned long int i;
  for(i=0; i<numOfNodes; i++) {
   	// update x,y coordinates
    nodeList[i].x += nodeList[i].dx;
    nodeList[i].y += nodeList[i].dy;
	
    // keep the nodes within the boundaries of the torus
    if (nodeList[i].x > TORUS_X) nodeList[i].x -= TORUS_X;
    if (nodeList[i].y > TORUS_Y) nodeList[i].y -= TORUS_Y;
    if (nodeList[i].x < 0) nodeList[i].x += TORUS_X;
    if (nodeList[i].y < 0) nodeList[i].y += TORUS_Y;
  }
}


/*
 * Cleanup function to be called before exit
 */
void cleanup(void) {
    double sum_Xi = 0;  // numerator of (1)
    double cp = 0;      // normalized packet cost
    double delta_i = 0; // stage duration
    double V_p = 0, C_p = 0, E_Xw = 0, E_Xb = 0, E_C = 0, E_D = 0;
    double Xw = 0, Xb = 0, Delta = 0;
    // allocate and initialize a histogram
	hist = gsl_histogram_alloc(HIST_BINS);
	gsl_histogram_set_ranges_uniform(hist, -M_PI, M_PI);
    // print the stages information into a file
    if (log_stages)
    {
        f_stages = fopen("log_stages.txt","w+");
        fprintf(f_stages,"lambda\tV0\tR0\ttw\tSimDur\tArea\tPotential\tBoundary\ta\tb\n%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.3f\t%c\t%c\t%.2f\t%.2f\t%u\n", 
                lambda, v0, r0, tw, simDur, (double)(TORUS_X*TORUS_Y), potential, boundary, a, b, stageOrTime);
        fprintf(f_stages,"i\tDelta_i (dur)\tTheta_i (dir)\tX_i (txDist)\tC_i (txCost)\tN_i (numHops)\n");
    }
    
	int i;
	for(i=0; i<numUsedStages; i++) {
        delta_i = stages[i].tEnd-stages[i].tStart;
        
        if (log_stages)
        {
            fprintf(f_stages,"%u\t%f\t%f\t%f\t%f\t%u\n",
					i,
					delta_i,
					stages[i].dir,
					stages[i].txDist,
					stages[i].txCost,
					stages[i].numHops);   
        }
        
        // calculate the performance metrics according to eqns. (1) and (2)
        sum_Xi += (v0 * delta_i * cos(stages[i].dir) + stages[i].txDist);
        cp += stages[i].txCost; // don't forget to divide by sum_Xi
        
        // intermediate metrics
        Xw += stages[i].txDist;
        Xb += v0 * delta_i * cos(stages[i].dir);
        Delta += delta_i;
        
        // update the histogram of Psi
		gsl_histogram_increment(hist, stages[i].dir);
	}
                
    if (log_stages)
    {
        fclose(f_stages);
    }

    // print the histogram
    f_psi = fopen("psi.txt","w+");
    //gsl_histogram_fprintf(f_psi, hist, "%g", "%g"); // segfault here!
	gsl_histogram_free(hist);
    fclose(f_psi);
    
    // metrics
    V_p = sum_Xi/Delta;
    C_p = cp/sum_Xi;
    
    // intermediate metrics
    E_Xw = Xw/numUsedStages;
    E_Xb = Xb/numUsedStages;
    E_C = cp/numUsedStages;
    E_D = Delta/numUsedStages;
    
    f_metrics = fopen("metrics.txt","w+");
    fprintf(f_metrics,"V_p Cp E_Xw E_Xb E_C E_D \n%f %f %f %f %f %f", 
        V_p, C_p, E_Xw, E_Xb, E_C, E_D);
    fclose(f_metrics);
    
    printf("\n-----------------------------------------------------------\n");
    printf("V_p Cp E_Xw E_Xb E_C E_D \n%f %f %f %f %f %f\n", 
        V_p, C_p, E_Xw, E_Xb, E_C, E_D);
    printf("-----------------------------------------------------------");
    printf("\nsum Delta = %f, numUsedStages = %u\n", Delta, numUsedStages);
	
	// releases the memory taken by the stages
	free(stages);
}

/*
 * Signal SIGINT (Ctrl+C) handler
 */
void ctrlC_handler(int signo) {
	// ignore the same signal while we are in its handler
	signal(signo, SIG_IGN);
    
    //discard the last stages since it's incomplete
    numUsedStages--;
	cleanup();
    
    
	printf("\nAddio...\n");
	exit(0);
}

/* 
 * Return the incidence angle of a node n that enters the forwarding region
 * of the carrier node c at location (x,y). The function acos(double x) 
 * returns in [0, pi] and x must be in [-1,1].
 * 
 * Remember that we are in a torus.
 */
double get_inc_angle(struct Node * n, struct Node * c)
{
	double xc = c->x;
	double yc = c->y;
	double xn = n->x;
	double yn = n->y;
	double dstnc = sqrt(dist_square(n,c));
	double incAng = NAN;
	
	if(fabs(xc-xn) > dstnc) {
		if(xc <= xn) 
			xn -= TORUS_X;
		else 
			xn += TORUS_X;
	}
	if(fabs(yc-yn) > dstnc) {
		if(yc <= yn) 
			yn -= TORUS_Y;
		else 
			yn += TORUS_Y;
	}
	xn -= xc;
	yn -= yc;
	
	// return the incidence angle in [-pi, pi]
	if ((xn >= 0 && yn >= 0) || (xn < 0 && yn >= 0)) 
		incAng = acos(xn/dstnc);	// 1st and 2nd quadrants
	else 
		incAng = -acos(xn/dstnc); 	// 3rd and 4th quadrants ((xn >= 0 && yn < 0)  || (xn < 0 && yn < 0 ))
	 
	return incAng;
}

/* 
 * Return a pointer to an unused stage_t structure
 */
stage_t *get_new_stage(stage_t *st, unsigned int *numAlloc, unsigned int *numUsed)
{
	stage_t *p = NULL;
	if(*numUsed == *numAlloc) {
		*numAlloc += NUM_ALLOC_STAGES;
		st = (stage_t*)realloc(st, (*numAlloc) * sizeof(stage_t));
		if (st == NULL) return p;
		
		// this is needed when realloc cannot just expand the current
		// block of memory, but it need to move everyting somewhere else.
		// When this happens, the pointer stages needs to be update, so
		// that it will point to the beginning of the newly allocated
		// block of memory.
		stages = st;
	}
	
	p = st + (*numUsed);
	(*numUsed)++;
	
	return p;
}


/*
 * Return 1 if the node pointed by n is in the FR of the node pointed by c.
 * Return 0, otherwise.
 */
int is_in_fr(struct Node *n, struct Node *c)
{
    /* 
     * If the FR is:
     *  - a disc of radius a: b(phi) = a, forall phi
     *  - a cardioid of radius 0.5: b(phi) = 1 + 0.5*cos(phi)
     *  - ellipse
     */ 
     
    double b_phi = 0;
    double phi = 0;
    
    if (boundary == 'd') {
        // disk
        b_phi = a;
    }
    else if (boundary == 'e') {
        // ellipse
        phi = get_inc_angle(n, c);
        b_phi = a*(1-b*b)/(1-b*cos(phi));
    }
    else if (boundary == 'c') {
        // cardioid
        phi = get_inc_angle(n, c);
        b_phi = a + b*cos(phi);
    }
    
    if (sqrt(dist_square(n, c)) <= b_phi)
        return 1;
    
    return 0;
}

/*
 * Return 1 if the node pointed by n1 is has a higher 
 * potential than the node pointed by n2.
 * 
 * Return 0, otherwise.
 */
int is_better(struct Node *n1, struct Node *n2)
{
    double u_n1 = 0;
    double u_n2 = 0;
    
    if (potential == 'a') {
        // Potential U(theta, r) = -|theta|
        u_n1 = -fabs(n1->d);
        u_n2 = -fabs(n2->d);   
    }
    // TODO else if ...
    
    return (u_n1 > u_n2);
}

/*
 * Return a new directin in [-pi, pi] distributed according to f_D(x)
 * tw is Theta_w (See paper)
 */
double get_new_direction(double tw)
{
    double d = tw/2;
    double x, X;
    if (tw >= 1.57 && tw < 1.58)    // pi/2
    {
        return ((double)rand()/RAND_MAX) * 2 * M_PI - M_PI;        
    }
    else
    {
        X = ((double)rand()/RAND_MAX) * 2 * M_PI - M_PI;
        x = fabs(X);
        while ((x>d && x<(M_PI_2-d))||(x>(M_PI_2+d) && x<(M_PI-d)))
        {
            X = ((double)rand()/RAND_MAX) * 2 * M_PI - M_PI;
            x = fabs(X);
        }
        return X;
    }
}

/*
 * Main function
 */
int main (int argc, char *argv[])
{
	// read command line arguments
	if (argc < 13) {
        printf("Usage: %s lambda v0 r0 tw simDur dT potential_func boudary_func a b stepOrTime seedRand\n", argv[0]);
        return 1;
    }
    lambda     = atof(argv[1]);
    v0         = atof(argv[2]);
    r0         = atof(argv[3]);
    tw         = atof(argv[4]);
    simDur     = atof(argv[5]);
    dT         = atof(argv[6]);
    potential  = *argv[7];
    boundary   = *argv[8];
    a          = atof(argv[9]);
    b          = atof(argv[10]);
    stageOrTime = atoi(argv[11]);    // 0 for time, 1 for stages
    log_stages = atoi(argv[12]);     // 1 for logging stages, 0 for not.
    seedRand   = atoi(argv[13]);     // -1 for using time(NULL)
   
    numOfNodes = (unsigned long int)(lambda * TORUS_X * TORUS_Y);
    dS = v0 * dT;
	
	// Local variables
	int progk = 1;                  // used to print the simulation progress
    struct Node *nodes;
    nodes = (struct Node *) malloc(numOfNodes * sizeof(struct Node));
	struct Packet pkt;              // the tagged packet
	
	signal(SIGINT, ctrlC_handler);	// signal handler
	
    // init random with a seed
    if (seedRand < 0) {
        srand(time(NULL));
    }
    else {
        srand((unsigned int)seedRand);
    }
		
	// Allocate stages. We allocate NUM_ALLOC_STAGES at a time to optimize the time spent
    // for memory allocation.
	stages = (stage_t*)malloc(NUM_ALLOC_STAGES * sizeof(stage_t));
	if (stages == NULL){
        printf("Cannot allocate memory!");
        exit(1);
    }
        
    numAllocStages = NUM_ALLOC_STAGES;
	numUsedStages = 0;
    
    // create a directory for the results
    char dirName[70];
    sprintf(dirName, "./lambda-%.3f_v0-%.2f_r0-%.2f_tw-%.2f_simDur-%.2f_pot-%c_bound-%c_a-%.2f_b-%.2f_stageOrTime-%u", 
            lambda, v0, r0, tw, simDur, potential, boundary, a, b, stageOrTime);
    mkdir(dirName, 0777);
    chdir(dirName);
         
    // print sim info
    printf("lambda %.3f, v0 %.2f, r0 %.2f, tw %.2f, simDur %.2f,\npot %c, bound %c, a %.2f, b %.2f, stageOrTime %u\n\n", 
            lambda, v0, r0, tw, simDur, potential, boundary, a, b, stageOrTime);
	
	//------------------------------------------------------------------
	//		Serious stuff starts here
	//------------------------------------------------------------------
	
	// Init nodes' position
	unsigned long int i;
	for(i=0; i<numOfNodes; i++)	{
		// generate nodes uniformly distributed on a torus 
		// x is random in [0.0, TORUS_W], inclusive
		// y is random in [0.0, TORUS_H], inclusive
		nodes[i].nodeId = i;
		nodes[i].x 	= ((double)rand()/RAND_MAX) * TORUS_X;			// initial x
		nodes[i].y 	= ((double)rand()/RAND_MAX) * TORUS_Y;			// initial y
		nodes[i].t 	= next_time(r0);							    // duration of the first leg
		nodes[i].d 	= get_new_direction(tw);	            // initial random direction in [-pi, pi]
		nodes[i].dx = dS * cos(nodes[i].d);					// x increment for the first leg
		nodes[i].dy = dS * sin(nodes[i].d);					// y increment for the first leg
	}
	
	// init packet
	pkt.carrierId = 0;
	pkt.d 		  = nodes[0].d;
	
	//------------------------------------------------------------------
	//						Simulation starts here
	//------------------------------------------------------------------
	double simTime = 0;
	double oldPktDir;
    
    // get a new stage and init its elements
	stage_t *currStage  = get_new_stage(stages, &numAllocStages, &numUsedStages);
	currStage->dir 		= pkt.d;
	currStage->tStart 	= simTime;
	currStage->numHops	= 0;
	currStage->txDist	= 0;
	currStage->txCost   = 0;
	
    // simulation loop
	while (1) {
		simTime+=dT;
		oldPktDir = pkt.d;
		
		// search for nodes that will have to change direction within the next dT
	    for (i=0; i<numOfNodes; i++) {
			nodes[i].t -= dT;
			if (nodes[i].t < 0) {   // it's time to change direction...
				// get a new waypoint
		        nodes[i].d 	= get_new_direction(tw);
		        nodes[i].t 	= next_time(r0);            // get an exp distributed r.v. sample
		        nodes[i].dx = dS * cos(nodes[i].d);
		        nodes[i].dy = dS * sin(nodes[i].d);
		        nodes[i].t -= dT;
			}
		}
		move_nodes((struct Node *)nodes); // move the nodes within the torus
		
	    // the packet has moved: update its direction (may not have changes)
		pkt.d = nodes[pkt.carrierId].d;
		
		double dstnc2 = 0;
	    
        // Here we do the routing. We get out of the while(1) when no more eligible nodes
        // are found. It can take forever.
	    while(1) {
			// to avoid an infinite loop (i.e., an eligible node is always found, this 
            // happens at large lambda or large forwarding region), use stageOrTime = 1
			
			// Find the best eligible in range of the carrier node
			int bestId = pkt.carrierId;
			for (i=0; i<numOfNodes; i++){
                
                // ------------ Optimization trick -------------------------------------//
                
                // The following block of code can be uncommented _ONLY_ if we are using a 
                // circular disc as forwarding region. It can be adapted to other shape of
                // FR.
                
                // skip nodes that are for sure not in range and the carrier
				// double absX = fabs(nodes[pkt.carrierId].x - nodes[i].x);
				// double absY = fabs(nodes[pkt.carrierId].y - nodes[i].y);
				// if (min(absX, TORUS_X-absX) > r || 
					// min(absY, TORUS_Y-absY) > r ||
					// pkt.carrierId == nodes[i].nodeId)
					// continue;
                    
                // ------------- End of optimization trick -----------------------------//

                // Is the node inside the FR of the carrier, AND its potential
                // is larger than the potential of the best node found so far (bestId)?
                if (is_in_fr(&nodes[i], &nodes[pkt.carrierId]) && is_better(&nodes[i], &nodes[bestId])) {
                    // we found a better eligible node
                    bestId = nodes[i].nodeId;
                }
			}
			
			// at this point the best node in range is found, if any.
			if (bestId != pkt.carrierId) {	// if there is a best node
				// update the tx progress, cost and num hops
				dstnc2 = dist_square(&nodes[pkt.carrierId], &nodes[bestId]);
				currStage->txDist += (sqrt(dstnc2) * cos(get_inc_angle(&nodes[bestId], &nodes[pkt.carrierId])));
				currStage->txCost += dstnc2;
				currStage->numHops++;
				
				// update packet
				pkt.carrierId = bestId;
				pkt.d = nodes[bestId].d;
			}
			else        // no eligible node is found...
				break;
		}
		
		//	condition to start a new stage
		if (oldPktDir != pkt.d)	{
			// end the current stage
			currStage->tEnd = simTime;
			
			// start a new one
			currStage  = get_new_stage(stages, &numAllocStages, &numUsedStages);
			if (currStage != NULL) {
				currStage->dir 		= pkt.d;
				currStage->tStart 	= simTime;
				currStage->txDist   = 0;
				currStage->txCost   = 0;
				currStage->numHops	= 0;
			}
		}
				
		// print simulation progress
		double p = (stageOrTime ? (numUsedStages/simDur * 100.0) : (simTime/simDur * 100.0));
		if(p >= 0.1*progk) {
			printf("\r%.1f%%",p);
			fflush(stdout);
			progk += 1;
		}
        
        // check if the simulation is over
        if (stageOrTime)
        {
            if (numUsedStages >= simDur)
                break;
        }
        else
        {
            if (simTime >= simDur)
                break;
        }
        
	}	// End of simulation
	
	// fix the last stage
	currStage -> tEnd = simTime;
	
    free(nodes);
	cleanup();
    
	return 0;
}
