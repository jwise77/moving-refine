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
 
int grid::AdjustVelocity(const float *velocity)
{
  
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  int i, j, k, dim, index, size;
  float TotalEnergy;
  
  for (dim = 0, size = 1; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* Find fields: density, total energy, velocity1-3. */

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }

  for (index = 0; index < size; index++) {
    TotalEnergy = 0.0;
    for (dim = 0; dim < GridRank; dim++) {
      BaryonField[Vel1Num+dim][index] -= velocity[dim];
      TotalEnergy += BaryonField[Vel1Num+dim][index] * BaryonField[Vel1Num+dim][index];
    }

    if (DualEnergyFormalism == TRUE)
      BaryonField[TENum][index] = 0.5*TotalEnergy + BaryonField[GENum][index];
    else
      BaryonField[TENum][index] = 0.5*TotalEnergy;
  } // ENDFOR index

  return SUCCESS;

}
