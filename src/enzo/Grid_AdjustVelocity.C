/***********************************************************************
/
/  GRID CLASS (CHANGE THE VELOCITY BY A SPECIFIED AMOUNT)
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
 
int grid::AdjustVelocity(float *velocity)
{

  int i, j, k, dim, index;
  float TotalEnergy;
  
  /* Find fields: density, total energy, velocity1-3. */

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      index = GRIDINDEX_NOGHOST(GridStartIndex[0], j, k);
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {

	TotalEnergy = 0.0;
	for (dim = 0; dim < GridRank; dim++) {
	  BaryonField[Vel1Num+dim][index] -= velocity[dim];
	  TotalEnergy += BaryonField[Vel1Num+dim][index] * BaryonField[Vel1Num+dim][index];
	}

	if (DualEnergyFormalism == TRUE)
	  BaryonField[TENum][index] = 0.5*TotalEnergy + BaryonField[GENum][index];
	else
	  BaryonField[TENum][index] = 0.5*TotalEnergy;
	  

      } // ENDFOR i
    } // ENDFOR j
  } // ENDFOR k

  return SUCCESS;

}
