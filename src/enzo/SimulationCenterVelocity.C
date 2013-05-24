/***********************************************************************
/
/  CENTER SIMULATION REFERENCE FRAME ON REFINE REGION
/
/  written by: John Wise
/  date:       May, 2013
/  modified1:  
/
/  PURPOSE:
/
************************************************************************/
#include "preincludes.h"
 
#ifdef USE_MPI
#include <mpi.h>
#endif
 
#include <stdio.h>

#include "EnzoTiming.h" 
#include "performance.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "CosmologyParameters.h"
#include "communication.h"
#include "CommunicationUtilities.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int SimulationCenterVelocity(LevelHierarchyEntry *LevelArray[], TopGridData *MetaData)
{

  int i, dim, level;
  LevelHierarchyEntry *ThisGrid;
  float velocity[MAX_DIMENSION];
  float total_mass = 0.0;

  if (HydroMethod == 2)
    ENZO_FAIL("CenterVelocityOnRefineRegion not tested with ZEUS.");

  /* Get units. */

  float LengthUnits, TimeUnits, TemperatureUnits, VelocityUnits, 
    DensityUnits;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, MetaData->Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  /* Determine mass-weighted velocity of the refine region */

  for (dim = 0; dim < MAX_DIMENSION; dim++)
    velocity[dim] = 0.0;

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
    for (ThisGrid = LevelArray[level]; ThisGrid; ThisGrid = ThisGrid->NextGridThisLevel) {

      ThisGrid->GridData->ComputeAverageVelocity(velocity, total_mass);

    } // ENDFOR grid
  } // ENDFOR level

  /* Divide out mass-weighting */

  for (dim = 0; dim < MAX_DIMENSION; dim++)
    velocity[dim] /= total_mass;

  if (debug)
    printf("Adjusting entire simulation by velocity (%g %g %g) cm/s\n",
	   velocity[0]*VelocityUnits, velocity[1]*VelocityUnits, velocity[2]*VelocityUnits);

  /* Adjust velocities in grids */

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
    for (ThisGrid = LevelArray[level]; ThisGrid; ThisGrid = ThisGrid->NextGridThisLevel) {

      ThisGrid->GridData->AdjustVelocity(velocity);

    } // ENDFOR grid
  } // ENDFOR level
  
  
  /* Turn off this flag because we don't want to recenter on a restart */

  CenterVelocityOnRefineRegion = INT_UNDEFINED;

  return SUCCESS;

}
