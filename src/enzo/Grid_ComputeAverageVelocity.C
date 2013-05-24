/***********************************************************************
/
/  GRID CLASS (COMPUTE THE AVERAGE VELOCITY)
/
/  written by: John Wise
/  date:       May, 2013
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/
************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
int grid::ComputeAverageVelocity(float *velocity, float &total_mass)
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  int i, j, k, dim, index, region;
  FLOAT CellVolume, CellMass;

  region = CenterVelocityOnRefineRegion;

  /* Check if the grid is in the specified refine region */

  for (dim = 0; dim < GridRank; dim++)
    if (GridLeftEdge[dim] < MultiRefineRegionLeftEdge[region][dim] ||
	GridRightEdge[dim] > MultiRefineRegionRightEdge[region][dim])
      return SUCCESS;

  /* Find fields: density, total energy, velocity1-3. */

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }

  CellVolume = 1.0;
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    CellVolume *= CellWidth[dim][0];

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      index = GRIDINDEX_NOGHOST(GridStartIndex[0], j, k);
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {

	CellMass = CellVolume * BaryonField[DensNum][index];
	for (dim = 0; dim < GridRank; dim++) {
	  velocity[dim] += BaryonField[Vel1Num+dim][index] * CellMass;
	  total_mass += CellMass;
	} // ENDFOR dim

      } // ENDFOR i
    } // ENDFOR j
  } // ENDFOR k

  return SUCCESS;

}
