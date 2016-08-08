
static char help[] = "Solves a tridiagonal linear system with KSP.\n\n";

/*T
   Concepts: KSP^solving a system of linear equations
   Processors: 1
T*/

/*
  Include "petscksp.h" so that we can use KSP solvers.  Note that this file
  automatically includes:
     petscsys.h       - base PETSc routines   petscvec.h - vectors
     petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners

  Note:  The corresponding parallel example is ex23.c
*/
#include <petscksp.h>
#include <stdlib.h>
#include <string.h>

char* concat(char *s1, char *s2)
{
    char *result = malloc(strlen(s1)+strlen(s2)+1);
    strcpy(result, s1);
    strcat(result, s2);
    return result;
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
    Vec             x;                      /* approx solution */
    Vec             b;                      /* RHS */
    Vec             u;                      /* exact solution */
    Mat             A;                      /* linear system matrix */
    KSP             ksp;                    /* linear solver context */
    PC              pc;                     /* preconditioner context */
    PetscReal       norm;                   /* norm of solution error */
    PetscRandom     rctx;                   /* random number generator context */
    PetscErrorCode  ierr;
    PetscInt        Si;
    PetscInt        its;                    /* number of iterations */
    PetscInt        R,C;                    /* index in matrix */
    PetscInt        wi,di;                  /* width of stencil along one direction, number of directions */
    PetscInt        Rstart,Rend;
    PetscInt        P[3]={0,0,0};           /* index in physical domain */
    PetscInt        D[3]={0,0,0};           /* displacement in physical domain */
    PetscScalar     neg_one = -1.0;
    PetscScalar     one = 1.0;
    PetscScalar     Avalue[4]={ 2.722222,-1.500000, 0.150000,-0.011111};
    PetscScalar     Lvalue[4]={ 1.983290,-0.729066, 0.073573,-0.005602};
    PetscBool       nonzeroguess = PETSC_FALSE;

    PetscInitialize(&argc,&args,(char*)0,help);

    //----------------------------------------------------------------------------------------------------------------------
    // Compute the matrix and right-hand-side vector that define the linear system, Ax = b.
    //----------------------------------------------------------------------------------------------------------------------

    for (Si=100; Si<101; Si++)
    {
        //----------------------------------------------------------------------------------------------------------------------
        // Set parameters.
        //----------------------------------------------------------------------------------------------------------------------
        PetscInt        N[3]={Si,Si,Si};     /* number of point along each direction */
        PetscInt        n=N[0]*N[1]*N[2];
        ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
        ierr = PetscOptionsGetBool(NULL,NULL,"-nonzero_guess",&nonzeroguess,NULL);CHKERRQ(ierr);
        //----------------------------------------------------------------------------------------------------------------------


        //----------------------------------------------------------------------------------------------------------------------
        // Create vectors.  Note that we form 1 vector from scratch and then duplicate as needed.
        //----------------------------------------------------------------------------------------------------------------------
        ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
        ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
        ierr = VecSetSizes(x,PETSC_DECIDE,n);CHKERRQ(ierr);
        ierr = VecSetFromOptions(x);CHKERRQ(ierr);
        ierr = VecDuplicate(x,&b);CHKERRQ(ierr);
        ierr = VecDuplicate(x,&u);CHKERRQ(ierr);         // if a vector exist
        //----------------------------------------------------------------------------------------------------------------------


        //----------------------------------------------------------------------------------------------------------------------
        // Create matrix.  When using MatCreate(), the matrix format can
        // be specified at runtime.
        // Performance tuning note:  For problems of substantial size,
        // preallocation of matrix memory is crucial for attaining good
        // performance. See the matrix chapter of the users manual for details.
        //----------------------------------------------------------------------------------------------------------------------
        ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
        ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
        ierr = MatSetFromOptions(A);CHKERRQ(ierr);
        ierr = MatMPIAIJSetPreallocation(A,19,NULL,0,NULL);CHKERRQ(ierr);
        ierr = MatSeqAIJSetPreallocation(A,19,NULL);CHKERRQ(ierr);
        ierr = MatSeqSBAIJSetPreallocation(A,1,19,NULL);CHKERRQ(ierr);
        ierr = MatSetUp(A);CHKERRQ(ierr);
        ierr = MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRQ(ierr);
        //----------------------------------------------------------------------------------------------------------------------


        //----------------------------------------------------------------------------------------------------------------------
        // Currently, all PETSc parallel matrix formats are partitioned by
        // contiguous chunks of rows across the processors.  Determine which
        // rows of the matrix are locally owned.
        //----------------------------------------------------------------------------------------------------------------------
        ierr = MatGetOwnershipRange(A,&Rstart,&Rend);CHKERRQ(ierr);
        //----------------------------------------------------------------------------------------------------------------------


        //----------------------------------------------------------------------------------------------------------------------
        // Create matrix for the 3-D, 19-points stencil in parallel
        //----------------------------------------------------------------------------------------------------------------------
        // Set matrix
        //----------------------------------------------------------------------------------------------------------------------
        for (R=Rstart; R<Rend; R++) 
        {
            P[0]=R/(N[1]*N[2]); P[1]=(R%(N[1]*N[2]))/N[2]; P[2]=R%N[2];
            for (di=0; di<3; di++)
            {
                // diagonal
                ierr=MatSetValues(A,1,&R,1,&R,&Avalue[0],ADD_VALUES);CHKERRQ(ierr);
                for (wi=1; wi<4; wi++)
                {
                    if (P[di]>=wi)         {D[0]=0; D[1]=0; D[2]=0; D[di]=-wi; C=(P[0]+D[0])*(N[1]*N[2])+(P[1]+D[1])*N[2]+(P[2]+D[2]); ierr=MatSetValues(A,1,&R,1,&C,&Avalue[wi],ADD_VALUES);CHKERRQ(ierr);}
                    if (P[di]<=N[di]-1-wi) {D[0]=0; D[1]=0; D[2]=0; D[di]= wi; C=(P[0]+D[0])*(N[1]*N[2])+(P[1]+D[1])*N[2]+(P[2]+D[2]); ierr=MatSetValues(A,1,&R,1,&C,&Avalue[wi],ADD_VALUES);CHKERRQ(ierr);}
                }
            }
        }
        //----------------------------------------------------------------------------------------------------------------------


        //----------------------------------------------------------------------------------------------------------------------
        // Assembly matrix
        //----------------------------------------------------------------------------------------------------------------------
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        //----------------------------------------------------------------------------------------------------------------------

        //----------------------------------------------------------------------------------------------------------------------
        // Check matrix
        //----------------------------------------------------------------------------------------------------------------------
        // MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        //----------------------------------------------------------------------------------------------------------------------


        //----------------------------------------------------------------------------------------------------------------------
        // Set exact solution; then compute right-hand-side vector.
        //----------------------------------------------------------------------------------------------------------------------
        ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx);CHKERRQ(ierr);
        ierr = PetscRandomSetFromOptions(rctx);CHKERRQ(ierr);
        ierr = VecSetRandom(u,rctx);CHKERRQ(ierr);
        ierr = PetscRandomDestroy(&rctx);CHKERRQ(ierr);
        ierr = MatMult(A,u,b);CHKERRQ(ierr);
        //----------------------------------------------------------------------------------------------------------------------

        //----------------------------------------------------------------------------------------------------------------------
        // Create the linear solver and set various options
        //----------------------------------------------------------------------------------------------------------------------

        //----------------------------------------------------------------------------------------------------------------------
        // Create linear solver context
        //----------------------------------------------------------------------------------------------------------------------
        ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
        //----------------------------------------------------------------------------------------------------------------------

        //----------------------------------------------------------------------------------------------------------------------
        // Set operators. Here the matrix that defines the linear system also serves as the preconditioning matrix.
        //----------------------------------------------------------------------------------------------------------------------
        ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
        ierr = KSPSetType(ksp,KSPCG);CHKERRQ(ierr);
        //----------------------------------------------------------------------------------------------------------------------

        //----------------------------------------------------------------------------------------------------------------------
        // Set linear solver defaults for this problem (optional).
        // - By extracting the KSP and PC contexts from the KSP context,
        //   we can then directly call any KSP and PC routines to set
        //   various options.
        // - The following four statements are optional; all of these
        //   parameters could alternatively be specified at runtime via
        //   KSPSetFromOptions();
        //----------------------------------------------------------------------------------------------------------------------
        ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
        ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
        ierr = KSPSetTolerances(ksp,1.e-6,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
        //----------------------------------------------------------------------------------------------------------------------

        //----------------------------------------------------------------------------------------------------------------------
        // Set runtime options, e.g.,
        //     -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
        // These options will override those specified above as long as
        // KSPSetFromOptions() is called _after_ any other customization
        // routines.
        //----------------------------------------------------------------------------------------------------------------------
        ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
        //----------------------------------------------------------------------------------------------------------------------

        if (nonzeroguess) 
        {
          PetscScalar p = 0.5;
          ierr = VecSet(x,p);CHKERRQ(ierr);
          ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);CHKERRQ(ierr);
        }

        //----------------------------------------------------------------------------------------------------------------------
        // Solve the linear system
        //----------------------------------------------------------------------------------------------------------------------
        ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
        //----------------------------------------------------------------------------------------------------------------------

        //----------------------------------------------------------------------------------------------------------------------
        // View solver info; we could instead use the option -ksp_view to
        // print this info to the screen at the conclusion of KSPSolve().
        //----------------------------------------------------------------------------------------------------------------------
        // ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        //----------------------------------------------------------------------------------------------------------------------

        //----------------------------------------------------------------------------------------------------------------------
        //                  Check solution and clean up
        //----------------------------------------------------------------------------------------------------------------------
        // Check the error
        //----------------------------------------------------------------------------------------------------------------------
        ierr = VecAXPY(x,neg_one,u);CHKERRQ(ierr);
        ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
        ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g, Iterations %D\n",(double)norm,its);CHKERRQ(ierr);
        //----------------------------------------------------------------------------------------------------------------------

        //----------------------------------------------------------------------------------------------------------------------
        // Free work space.  All PETSc objects should be destroyed when they are no longer needed.
        //----------------------------------------------------------------------------------------------------------------------
        ierr = VecDestroy(&x);   CHKERRQ(ierr); 
        ierr = VecDestroy(&u);   CHKERRQ(ierr);
        ierr = VecDestroy(&b);   CHKERRQ(ierr); 
        ierr = MatDestroy(&A);   CHKERRQ(ierr);
        ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
        //----------------------------------------------------------------------------------------------------------------------
    }

    //----------------------------------------------------------------------------------------------------------------------
    // Always call PetscFinalize() before exiting a program.  This routine
    //   - finalizes the PETSc libraries as well as MPI
    //   - provides summary and diagnostic information if certain runtime
    //     options are chosen (e.g., -log_summary).
    //----------------------------------------------------------------------------------------------------------------------
    ierr = PetscFinalize();
    return 0;
}
