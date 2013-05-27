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

inline bool inside_multi_refine_region(const FLOAT x, const FLOAT y, const FLOAT z, 
				       const int region)
{
  bool inside = true;
  inside &= (x >= MultiRefineRegionLeftEdge[region][0] &&
	     x <= MultiRefineRegionRightEdge[region][0]);
  if (inside) {
    inside &= (y >= MultiRefineRegionLeftEdge[region][1] &&
	       y <= MultiRefineRegionRightEdge[region][1]);
    if (inside)
      inside &= (z >= MultiRefineRegionLeftEdge[region][2] &&
		 z <= MultiRefineRegionRightEdge[region][2]);
  }
  return inside;
}
 
int grid::ComputeAverageVelocity(float *velocity, float &total_mass)
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  int i, j, k, dim, index, region;
  FLOAT x, y, z, CellVolume, CellMass;

  region = CenterVelocityOnRefineRegion;

  /* Check if the grid is in the specified refine region */

  bool contained = true;
  for (dim = 0; dim < GridRank; dim++)
    contained &= (GridLeftEdge[dim] < MultiRefineRegionLeftEdge[region][dim] &&
		  GridRightEdge[dim] > MultiRefineRegionRightEdge[region][dim]);
  bool outside = false;
  for (dim = 0; dim < GridRank; dim++)
    outside |= (GridRightEdge[dim] < MultiRefineRegionLeftEdge[region][dim] ||
		GridLeftEdge[dim] > MultiRefineRegionRightEdge[region][dim]);
  
  if (!contained && outside)
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
    z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
      index = GRIDINDEX_NOGHOST(GridStartIndex[0], j, k);
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {

	if (BaryonField[NumberOfBaryonFields][index] < 1.0) {
	  x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	  if (inside_multi_refine_region(x, y, z, region)) {
	    CellMass = CellVolume * BaryonField[DensNum][index];
	    for (dim = 0; dim < GridRank; dim++) {
	      velocity[dim] += BaryonField[Vel1Num+dim][index] * CellMass;
	    } // ENDFOR dim
	    total_mass += CellMass;
	  } // ENDIF inside
	} // ENDIF no subgrid

      } // ENDFOR i
    } // ENDFOR j
  } // ENDFOR k

  return SUCCESS;

}
