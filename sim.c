/* Simulator.
 * 
 * Nodes moves on a torus.
 * 
 * For details see routing.pdf.
 * The program writes a log file that contains the histogram of the 
 * direction of the packet: each line corresponds to one bin of the 
 * histogram, in this way we avoid to write a huge file.
 * 
 * See http://www.gnu.org/software/gsl/manual/html_node/Histograms.html
 * for details.
 * 
 * Build and run with: 
 * 		gcc -Wall -lgsl -lgslcblas -lm -o sim sim.c && ./sim
 * 
 * How to read the histogram file in Matlab:
 * fp = fopen('log.txt');
 * a = textscan(fp,'%f %f %f'); % a is a cell 
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
#define TORUS_X			  	200
#define TORUS_Y			    200
#define HIST_BINS			101 // odd number so 0 is included
#define ALLOC_MIN_SIZE  	10
#define NUM_ALLOC_STAGES  	100
#define min(a, b) (((a) < (b)) ? (a) : (b))
//#define DEBUG

// Typedefs
struct Node {
	int nodeId;	// id of the node
	double x;		// current x-position
	double y;		// current y-position
	double t;		// time to travel before next turn
	double d;		// direction of travel
	double dx;		// x increment during a leg
	double dy;		// y increment during a leg
};

struct Packet {
	double d;		// current direction of travel
	double x;		// current x-position
	double y;		// current y-position
	int carrierId; 		// id of the current carrier node
	size_t numBanned;	// size of the banned list
	int * banned; 		// pointer to an array of banned nodes
};

typedef struct {
	double dir;		 // direction
	double tStart;	 // start time
	double tEnd;	 // end time
	double txDist;	 // distance covered by tx (can be zero)
	double txCost;   // tx cost (can be zero)
	unsigned int numHops;	// number of hops
}stage_t;

// Global varialbles
unsigned int numOfNodes;
double r;
double r2;
double v0;
double r0;
double simDur;
double dT;
double dS;
int rule;
double tau;
FILE * fp;
FILE * fp3;
gsl_histogram * hist = NULL;
stage_t *stages = NULL;
unsigned int numAllocStages;	// sizeof(unsigned int) = 4
unsigned int numUsedStages;	
	

/*
 * Generate samples of an exponentially distribute random variable
 * with a given rate.
 */
double next_time(double rate)
{	
	double r = 0;
	while (r == 0) {	
		r = (double)rand()/RAND_MAX;
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
  int i;
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
 * Return 1 if array arr contains element val,
 * return 0 otherwise.
 */
int is_in_array(int val, int *arr, size_t len)
{
	if(arr == NULL) return 0;
	
	int i;
	for(i=0; i<len; i++) {
		if(arr[i] == val) return 1;
	}
	return 0;
}

/* 
 *  Create a list of nodes banned by the routing rule.
 *  A node is banned if it is in the forw. reg. of the carrier.
 * 	In this way we avoid multihop.
 */
void create_ban_list(struct Packet * pkt, struct Node * nodeList)
{
	// release the previous array and allocate a new one
	free(pkt->banned);
	pkt->numBanned = 0;
	pkt->banned = malloc(ALLOC_MIN_SIZE * sizeof(int));
	
	int i;
	for(i=0; i<numOfNodes; i++)	{
		// To be banned a node has to:
		// 1) not be the current carrier
		// 2) being in the tx range of the current carrier
		
		if( nodeList[i].nodeId != pkt->carrierId &&
			dist_square(&nodeList[i], &nodeList[pkt->carrierId]) < r2) {
			// add node to the banned list
			pkt->numBanned += 1;
			pkt->banned[pkt->numBanned - 1] = nodeList[i].nodeId;
			
			if(pkt->numBanned % ALLOC_MIN_SIZE == 0) {
				pkt->banned = realloc(pkt->banned, (pkt->numBanned + ALLOC_MIN_SIZE)*sizeof(int));
			}
		}
	}
}

/*
 * Cleanup function to be called before exit
 */
void cleanup(void) {
	
	// write and delete the histogram
	gsl_histogram_fprintf(fp, hist, "%g", "%g");
	gsl_histogram_free(hist);
  
	// close the files
	fclose(fp);
	
	// print the stages information on a file
	printf("used = %u\talloc = %u\n", numUsedStages, numAllocStages);
	fprintf(fp3,"num\tdur\tdir\ttxDist\ttxCost\tnumHops\n");
	int i;
	for(i=0; i<numUsedStages; i++) {
		fprintf(fp3,"%u\t%f\t%f\t%f\t%f\t%u\n",
					i,
					stages[i].tEnd-stages[i].tStart,
					stages[i].dir,
					stages[i].txDist,
					stages[i].txCost,
					stages[i].numHops);
	}
	fclose(fp3);
	
	// releases the memory taken by the stages
	free(stages);
}

/*
 * Signal SIGINT (Ctrl+C) handler
 */
void ctrlC_handler(int signo) {
	// ignore the same signal while we are in its handler
	signal(signo, SIG_IGN);		
	printf("\nAddio...\n");
	cleanup();
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
 * Implementation of the routing rules.
 * Return 1 if the encountered node is eligible, 0 otherwise
 */
int routingRule(int rule, int id, int bestId, int carrId, struct Node * nList){
	
	int ret = 0;
	
	// Routing 1
	if (rule == 1 && 
		fabs(nList[id].d) < fabs(nList[bestId].d))
		ret = 1;
	
	// Routing 2
	else if (rule == 2 && 
			 fabs(get_inc_angle(&nList[id], &nList[carrId])) <= M_PI_2 &&
			 fabs(nList[id].d) < fabs(nList[bestId].d))
		ret = 1;

	// Routing 3
	else if (rule == 3 && 
			 fabs(get_inc_angle(&nList[id], &nList[carrId])) <= M_PI_2 &&
			 dist_square(&nList[id], &nList[carrId]) > dist_square(&nList[bestId], &nList[carrId]))
		ret = 1;
	
	// Routing 4
	else if (rule >= 4 &&
			fabs(nList[id].d) <= tau &&
			fabs(nList[id].d) < fabs(nList[bestId].d))
		ret = 1;
	
	else
		ret = 0;
	
	return ret;
}

/*
 * Main function
 */
int main (int argc, char *argv[])
{
	// read command line arguments
	if (argc < 8) {
            printf("Usage: %s numOfNodes r v0 r0 simDur dT RULE_NUM\n",argv[0]);
            return 1;
    }
    numOfNodes = atoi(argv[1]);
    r 		   = atof(argv[2]);
    v0         = atof(argv[3]);
    r0         = atof(argv[4]);
    simDur     = atof(argv[5]);
    dT         = atof(argv[6]);
    rule 	   = atoi(argv[7]);
    dS 		   = v0 * dT;
    if(rule >= 4) {
		if(argc == 9)
			tau = atof(argv[8]);
		else {	
			printf("You forgot tau\n");
			exit(1);
		}
	}
	r2 = r*r;
	
	// Local variables
	int progk = 1;
	struct Node nodes[numOfNodes];
	struct Packet pkt;
	
	signal(SIGINT, ctrlC_handler);	// signal handler
	srand(time(NULL));				// init rand with a seed
	//~ srand(1);

	// allocate and initialize a histogram
	hist = gsl_histogram_alloc(HIST_BINS);
	gsl_histogram_set_ranges_uniform(hist, -M_PI, M_PI);
		
	// stages
	stages = (stage_t*)malloc(NUM_ALLOC_STAGES * sizeof(stage_t));
	numAllocStages = NUM_ALLOC_STAGES;	// sizeof(unsigned int) = 4
	numUsedStages = 0;
    
    // create a directory for the results
    char dirName[70];
    sprintf(dirName, "./NODES-%u_R-%.2f_V0-%.2f_R0-%.2f_SIMDUR-%.2f_RULE-%u", numOfNodes, r, v0, r0, simDur, rule);
    mkdir(dirName, 0777);
    chdir(dirName);
    
    fp  = fopen("log.txt","w+");
	fp3 = fopen("log_stages.txt","w+");
	
	//------------------------------------------------------------------
	//		Serious stuff starts here
	//------------------------------------------------------------------
	
	// init the direction array: nodes' directions are discretized
	unsigned int numDir = 1000;	// don't change this
	double drctn[numDir+1];
	double dstep = 2*M_PI/numDir;
	int dirSize = sizeof(drctn)/sizeof(drctn[0]);
	int i;
	for (i=0; i<=numDir; i++){
		drctn[i] = -M_PI + i*dstep;
	}
	
	// Init nodes' position
	//~ int i;
	for(i=0; i<numOfNodes; i++)	{
		// generate nodes uniformly distributed on a torus 
		// x is random in [0.0, TORUS_W], inclusive
		// y is random in [0.0, TORUS_H], inclusive
		nodes[i].nodeId = i;
		nodes[i].x 	= ((double)rand()/RAND_MAX) * TORUS_X;			// initial x
		nodes[i].y 	= ((double)rand()/RAND_MAX) * TORUS_Y;			// initial y
		nodes[i].t 	= next_time(r0);									// duration of the first leg
		// nodes[i].d 	= ((double)rand()/RAND_MAX) * 2 * M_PI - M_PI;	// initial direction
		nodes[i].d 	= drctn[rand() % dirSize];	// initial direction
		nodes[i].dx = dS * cos(nodes[i].d);					// x increment for the first leg
		nodes[i].dy = dS * sin(nodes[i].d);					// y increment for the first leg
	}
	
	// init packet
	pkt.carrierId 	= 0;
	pkt.d 			= nodes[0].d;
	pkt.numBanned 	= 0;
	pkt.banned 		= malloc(pkt.numBanned*sizeof(pkt.carrierId));
	
	// discard nodes initially located inside the forwarding region
	//create_ban_list(&pkt, nodes);
	
	//------------------------------------------------------------------
	//						Simulation starts here
	//------------------------------------------------------------------
	double simTime = 0;
	double oldPktDir;
	stage_t *currStage  = get_new_stage(stages, &numAllocStages, &numUsedStages);
	currStage->dir 		= pkt.d;
	currStage->tStart 	= simTime;
	currStage->numHops	= 0;
	currStage->txDist	= 0;
	currStage->txCost	= 0;
	
	while (simTime<simDur) {
		simTime+=dT;
		oldPktDir = pkt.d;
		
		// move the nodes
	    for (i=0; i<numOfNodes; i++) {	
			nodes[i].t -= dT;
			if (nodes[i].t < 0) {
				// get a new waypoint
		        // nodes[i].d 	= ((double)rand()/RAND_MAX) * 2 * M_PI - M_PI;
		        nodes[i].d 	= drctn[rand() % dirSize];
		        nodes[i].t 	= next_time(r0);
		        nodes[i].dx = dS * cosf(nodes[i].d);
		        nodes[i].dy = dS * sinf(nodes[i].d);
		        nodes[i].t -= dT;
			}
		}
		move_nodes((struct Node *)nodes); // move the nodes
		
	    // the packet has moved: update its direction
		pkt.d = nodes[pkt.carrierId].d;
		
		double dstnc2 = 0;
	    
	    while(1) {
			// if the packet is in the best direction there is no way
			// to transmit it to anybody else.
			if (pkt.d == drctn[numDir/2])
				break;
			
			// Find the best node range of the carrier node
			int bestId = pkt.carrierId;
			for (i=0; i<numOfNodes; i++){
				// skip nodes that are for sure not in range and the carrier
				double absX = fabs(nodes[pkt.carrierId].x - nodes[i].x);
				double absY = fabs(nodes[pkt.carrierId].y - nodes[i].y);
				if (min(absX, TORUS_X-absX) > r || 
					min(absY, TORUS_Y-absY) > r ||
					pkt.carrierId == nodes[i].nodeId)
					continue;

				dstnc2 = dist_square(&nodes[pkt.carrierId], &nodes[i]);
				if (dstnc2 <= r2 && routingRule(rule, i, bestId, pkt.carrierId, nodes)) {
					#ifdef DEBUG
						printf("t = %f, better node %d\n", simTime, nodes[i].nodeId);			
					#endif //DEBUG
					bestId = nodes[i].nodeId;
				}
				
				#ifdef DEBUG
					if (dstnc2 <= r2)
						printf("t = %f, node %d in range, xn = %f, xc = %f, yn = %f, yc = %f, incAng = %f, dist = %f, dir = %f, mydir = %f\n",
								simTime, nodes[i].nodeId, nodes[i].x, nodes[pkt.carrierId].x,
								nodes[i].y, nodes[pkt.carrierId].y,
								get_inc_angle(&nodes[i],&nodes[pkt.carrierId]),
								sqrt(dstnc2), nodes[i].d, nodes[pkt.carrierId].d);
				#endif //DEBUG
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
				
				#ifdef DEBUG
					printf("t = %f, best node: %d\n-------------------------\n",
							simTime, nodes[bestId].nodeId);
					printf("Press Any Key to Continue\n");  
					getchar();  
				#endif //DEBUG
			}
			else
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
			
	    // update the histogram: it must be done once per timestep
		gsl_histogram_increment(hist, pkt.d);
			
		// print simulation progress
		double p = simTime/simDur * 100.0;
		if(p >= 0.1*progk) {
			printf("\r%.1f%%",p);
			fflush(stdout);
			progk += 1;
		}
	}	// End of simulation
	
	// fix the last stage
	currStage -> tEnd = simTime;
	
	cleanup();
	return 0;
}
