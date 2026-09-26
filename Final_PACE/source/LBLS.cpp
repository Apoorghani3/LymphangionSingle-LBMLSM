// LBLS.cpp : Defines the entry point for the console application.
// (c) Alexander Alexeev, 2006 
//

#include "clVector.h"
#include "StdAfx.h"
#include "clConstants.h"
#include "clMethod.h"
#include "clVariables.h"
#include "clProblem.h"
#include "clOutput.h"
#include "clDisplay.h"
clConstants c;
clVariables v;
clVariablesLS vls;
clVariablesLB vlb;
clVariablesTr vtr;
clMethod m;
clProblem p;


int main(int argc, char* argv[])
{

#ifdef USE_OPENMP
#ifdef NUM_THREADS
	omp_set_num_threads(NUM_THREADS);
#endif //NUM_THREADS
#pragma omp parallel default(shared)
#pragma omp master
	printf("OpenMP threads: %d\n", omp_get_num_threads());
#endif //USE_OPENMP

	p.ini();
	p.run();
	return 0;
}
