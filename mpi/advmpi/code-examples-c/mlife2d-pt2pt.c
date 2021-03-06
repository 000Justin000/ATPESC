/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  (C) 2004 by University of Chicago.
 *      See COPYRIGHT in top-level directory.
 */

#include <mpi.h>

#include "mlife2d.h"

int MLIFE_ExchangeInitPt2pt( MLIFEPatchDesc *patch, int **m1, int **m2, void *privateData )
{
    *(void **)privateData = 0;
    return 0;
}

int MLIFE_ExchangeEndPt2pt( void *privateData )
{
    return 0;
}

int MLIFE_ExchangePt2pt( MLIFEPatchDesc *patch, int **matrix, MLIFETiming *timedata, void *privateData )
{
    MPI_Request reqs[4];
    MPI_Comm    comm = patch->comm;
    static MPI_Datatype type = MPI_DATATYPE_NULL;
    int LRows = patch->lni;
    int LCols = patch->lnj;

    /* Send and receive boundary information */
    /* For a non-contigous data, this is useful */
    if (type == MPI_DATATYPE_NULL) 
    {
        MPI_Type_vector(LRows, 1, LCols+2, MPI_INT, &type);
        MPI_Type_commit(&type);
    }

    /* first, move the left, right edges */
    MPI_Irecv(&matrix[1][0], 	   1, type, patch->left,  0, comm, reqs+0);
    MPI_Irecv(&matrix[1][LCols+1], 1, type, patch->right, 0, comm, reqs+1);
    MPI_Isend(&matrix[1][1], 	   1, type, patch->left,  0, comm, reqs+2);
    MPI_Isend(&matrix[1][LCols],   1, type, patch->right, 0, comm, reqs+3);

    /* We need to wait on these for the trick that we use to move the diagonal terms to work */
    MPI_Waitall( 4, reqs, MPI_STATUSES_IGNORE );

    /* move the top, bottom edges (including diagonals) */
    MPI_Irecv(&matrix[0][0],       LCols+2, MPI_INT, patch->up,   0, comm, reqs+0);
    MPI_Irecv(&matrix[LRows+1][0], LCols+2, MPI_INT, patch->down, 0, comm, reqs+1);
    MPI_Isend(&matrix[1][0],       LCols+2, MPI_INT, patch->up,   0, comm, reqs+2);
    MPI_Isend(&matrix[LRows][0],   LCols+2, MPI_INT, patch->down, 0, comm, reqs+3);

    MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);

    return MPI_SUCCESS;
}
