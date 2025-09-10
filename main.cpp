#include <TPZYSMPMatrix.h>
#include <pzstepsolver.h>
#include <fstream>
#include <iostream>
#include <pzgmesh.h> //for TPZGeoMesh
#include "TPZVTKGeoMesh.h"
#include "TPZAnalyticSolution.h"
#include "DarcyFlow/TPZHybridDarcyFlow.h"//I can't invoke only the TPZDarcyFlow class
#include "TPZLinearAnalysis.h"
#include "TPZSSpStructMatrix.h"
#include "pzgeoquad.h"
#include "tpzgeoelrefpattern.h"

TPZGeoMesh* CreateQuadLShapeMesh(TPZVec<int>& bcids);
TPZGeoMesh* CreateTriangLShapeMesh(TPZVec<int>& bcids);
void UniformRefinement(int nDiv, TPZGeoMesh* gmesh);
TPZCompMesh* CreateCMeshH1(TPZGeoMesh* geomesh,TLaplaceExample1* exactsol, int intorder, int porder);


int main(){
    
    int porder = 1;
    int nref = 4; // Number of refinements to be applied to the initial mesh
    int nthreads = 0;
    int integrationorder = 11;
    std::string topology = "Quadrilateral"; //Triangular, Quadrilateral

    std::string problemname = "ESinSin";//ESinSin,ESinMark,EConst,EBubble2D,ESteepWave;
    TLaplaceExample1 aux, exact;
    exact.fExact = aux.ESinSin;//ESinMark//ESinSin//ESinSinDirNonHom

    // Create geometric mesh
    TPZManVector<int, 8> Lshape_bcids(8, -1);
    TPZGeoMesh *gmesh = nullptr;
    
    if (topology == "Quadrilateral"){
    gmesh = CreateQuadLShapeMesh(Lshape_bcids);
    //gmesh->Print();
    }
    
    if (topology == "Triangular"){
    gmesh = CreateTriangLShapeMesh(Lshape_bcids);
    //gmesh->Print();
    }
    
    std::ofstream salida("mallageometrica.txt");
    gmesh->Print(salida);
    
    std::ofstream out("mallageom.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    
    UniformRefinement(nref, gmesh);
    std::ofstream out2("mallageomrefinada.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2);
    
    //Create computational mesh
    TPZCompMesh* cmesh = nullptr;
    cmesh = CreateCMeshH1(gmesh, &exact, integrationorder, porder);
    std::ofstream out3("mallacomputac.txt");
    cmesh->Print(out3);
    
    //Assembly and resolution of the linear system
    TPZLinearAnalysis an(cmesh);
    
    #ifdef PZ_USING_MKL
        TPZSSpStructMatrix<STATE> strmat(cmesh);
        strmat.SetNumThreads(0);
        //strmat.SetDecomposeType(ELDLt);
    #else
    //  TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(cmeshH1);
    //  strmat.SetNumThreads(0);
        TPZSkylineStructMatrix<STATE> strmat(cmeshH1);
        strmat.SetNumThreads(0);
    #endif
    
    std::set<int> matids={1,-1};
    strmat.SetMaterialIds(matids);
    an.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;
    an.Assemble();
    an.Solve();

    int64_t nelem = cmesh->NElements();
    cmesh->LoadSolution(cmesh->Solution());
    cmesh->ExpandSolution();
    cmesh->ElementSolution().Redim(nelem, 10);
    
    TPZManVector<REAL,3> error;
    std::ofstream anPostProcessFile("postprocess.txt");
    an.PostProcess(error,anPostProcessFile);
    
    //cmesh->Solution().Print("Solution");
    
    std::cout << "NDofs: " << cmesh->NEquations() << '\n';
    std::cout << "\nApproximation error:\n";
    std::cout << "H1 Norm = " << error[0]<<'\n';
    std::cout << "L2 Norm = " << error[1]<<'\n';
    std::cout << "H1 Seminorm = " << error[2] << "\n\n";
    
    //Log approximation errors in a file
    std::ofstream fileouput;

    fileouput.open("TrueErrors.txt",std::ios::app);
    fileouput << std::setw(15) <<"Problem name" << std::setw(15)<<"p-order" << std::setw(15) <<"DOF's" <<std::setw(15) <<"H1-error" << std::setw(15)<< "L2-error" <<std::setw(15) << "L2-seminorm" << std::endl;

    fileouput << std::setw(15) << problemname;
    fileouput << std::setw(15) << porder;
    fileouput << std::setw(15) << cmesh->NEquations();
    fileouput << std::setw(15) << error[0]; // H1-norm
    fileouput << std::setw(15) << error[1]; // L2-norm
    fileouput << std::setw(15) << error[2] << std::endl; // L2-seminorm
    fileouput.close();
    
    
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Solution");
    vecnames.Push("Derivative");
    vecnames.Push("Flux");
    scalnames.Push("ExactSolution");
    vecnames.Push("ExactFlux");
    int dim = cmesh->Reference()->Dimension();
    
    std::string plotname="poissonSolution.vtk";
    int resolution = 3;
    an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
    an.PostProcess(resolution,dim);
    
    return 0;
}


TPZGeoMesh* CreateQuadLShapeMesh(TPZVec<int>& bcids) {

    TPZGeoMesh* gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    int matID = 1;

    // Creates matrix with node coordinates
    const int NodeNumber = 8;
    REAL coordinates[NodeNumber][3] = {
            {0.,  0., 0.},
            {1.,  0., 0.},
            {1.,  1., 0.},
            {0.,  1., 0.},
            {-1., 1., 0.},
            {-1., 0., 0.},
            {-1.,-1., 0.},
            {0., -1., 0.}
    };

    // Inserts coordinates in the TPZGeoMesh object
    for (int i = 0; i < NodeNumber; i++) {
        int64_t nodeID = gmesh->NodeVec().AllocateNewElement();

        TPZVec<REAL> nodeCoord(3);
        nodeCoord[0] = coordinates[i][0];
        nodeCoord[1] = coordinates[i][1];
        nodeCoord[2] = coordinates[i][2];

        gmesh->NodeVec()[nodeID] = TPZGeoNode(i, nodeCoord, *gmesh);
    }

    // Creates 2D elements
    TPZManVector<int64_t> nodeIDs(4);
    for (int i = 0; i < 3; i++) {
        nodeIDs[0] = 0;
        nodeIDs[1] = (2 * i + 1) % NodeNumber;
        nodeIDs[2] = (2 * i + 2) % NodeNumber;
        nodeIDs[3] = (2 * i + 3) % NodeNumber;
        new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodeIDs, matID, *gmesh);
    }

    // Creates line elements where boundary conditions will be inserted
    nodeIDs.Resize(2);
    for (int i = 0; i < NodeNumber; i++) {
        nodeIDs[0] = i % NodeNumber;
        nodeIDs[1] = (i + 1) % NodeNumber;
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodeIDs, bcids[i], *gmesh);
    }

    gmesh->BuildConnectivity();

    return gmesh;

}

TPZGeoMesh* CreateTriangLShapeMesh(TPZVec<int>& bcids){
    
    TPZGeoMesh* gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    int matID = 1;
    
    // Creates matrix with quadrilateral node coordinates.
    const int NodeNumber = 8;
    REAL coordinates[NodeNumber][3] = {
        {0., 0., 0.},
        {1., 0., 0.},
        {1., 1., 0.},
        {0., 1., 0.},
        {-1.,1.,0.},
        {-1.,0.,0.},
        {-1.,-1.,0.},
        {0.,-1.,0.}
    };
    
    // Inserts coordinates in the TPZGeoMesh object.
    for(int i = 0; i < NodeNumber; i++) {
        int64_t nodeID = gmesh->NodeVec().AllocateNewElement();
        
        TPZVec<REAL> nodeCoord(3);
        nodeCoord[0] = coordinates[i][0];
        nodeCoord[1] = coordinates[i][1];
        nodeCoord[2] = coordinates[i][2];
        
        gmesh->NodeVec()[nodeID] = TPZGeoNode(i, nodeCoord, *gmesh);
    }
    
    // Creates triangular elements.
    int64_t index =0;
    TPZManVector<int64_t> nodeIDs(3);
    //El 0
    nodeIDs[0] = 0;
    nodeIDs[1] = 1;
    nodeIDs[2] = 3;
    gmesh->CreateGeoElement(ETriangle, nodeIDs, matID, index);
    index++;
    
    //El 1
    nodeIDs[0] = 2;
    nodeIDs[1] = 3;
    nodeIDs[2] = 1;
    gmesh->CreateGeoElement(ETriangle, nodeIDs, matID, index);
    index++;
    //El 2
    nodeIDs[0] = 3;
    nodeIDs[1] = 4;
    nodeIDs[2] = 0;
    gmesh->CreateGeoElement(ETriangle, nodeIDs, matID, index);
    index++;
    //El 3
    nodeIDs[0] = 5;
    nodeIDs[1] = 0;
    nodeIDs[2] = 4;
    gmesh->CreateGeoElement(ETriangle, nodeIDs, matID, index);
    index++;
    
    //El 4
    nodeIDs[0] = 0;
    nodeIDs[1] = 5;
    nodeIDs[2] = 7;
    gmesh->CreateGeoElement(ETriangle, nodeIDs, matID, index);
    index++;
    //El 6
    nodeIDs[0] = 6;
    nodeIDs[1] = 7;
    nodeIDs[2] = 5;
    gmesh->CreateGeoElement(ETriangle, nodeIDs, matID, index);
    index++;
    
    // Creates line elements where boundary conditions will be inserted.
    nodeIDs.Resize(2);
    
    for (int i = 0; i < NodeNumber-1; i++) {
        
        nodeIDs[0] = i;
        
        nodeIDs[1] = (i + 1);
        std::cout<<"xo "<<nodeIDs[0]<<" x1 "<<nodeIDs[1]<<" bcid "<<bcids[i]<< "\n";
        
        gmesh->CreateGeoElement(EOned, nodeIDs, bcids[i], index);
    }
    index ++;
    
    nodeIDs[0] = 7;
    nodeIDs[1] = 0;
    std::cout<<"xo "<<nodeIDs[0]<<" x1 "<<nodeIDs[1]<<" bcid "<<bcids[NodeNumber-1]<< "\n";
    
    gmesh->CreateGeoElement(EOned, nodeIDs, bcids[NodeNumber-1], index);
    
    
    gmesh->BuildConnectivity();
    
    return gmesh;
    
}

void UniformRefinement(int nDiv, TPZGeoMesh* gmesh) {
    
    TPZManVector<TPZGeoEl*> children;
    for (int division = 0; division < nDiv; division++) {
        
        int64_t nels = gmesh->NElements();
        
        for (int64_t elem = 0; elem < nels; elem++) {
            
            TPZGeoEl* gel = gmesh->ElementVec()[elem];
            
            if (!gel || gel->HasSubElement()) continue;
            if (gel->Dimension() == 0) continue;
            gel->Divide(children);
        }
    }
}


TPZCompMesh* CreateCMeshH1(TPZGeoMesh* geomesh,TLaplaceExample1* exactsol, int intorder, int porder) {

    TPZCompMesh* cmesh = new TPZCompMesh(geomesh);
    TPZDarcyFlow* mat = 0;
    int dirichlet = 0;
    int neumann = 1;
    
    int matid = 1;
    int bcmatid = -1;
    
    int dim = geomesh->Dimension();

    TPZDarcyFlow *mix = new TPZDarcyFlow(matid, cmesh->Dimension());
    mix->SetExactSol(exactsol->ExactSolution(), intorder);
    mix->SetForcingFunction(exactsol->ForceFunc(), intorder);

    if (!mat) mat = mix;
    cmesh->InsertMaterialObject(mix);
        
    TPZFNMatrix<1, REAL> val1(1, 1, 0.);
    TPZManVector<STATE, 2> val2(1, 0.); //Dirichlet
    int bctype = 0; //Dirichlet
    auto *bc = mat->CreateBC(mat, bcmatid, bctype, val1, val2);
    bc->SetForcingFunctionBC(exactsol->ExactSolution(),intorder);
    cmesh->InsertMaterialObject(bc);
        

    cmesh->SetDefaultOrder(porder);//ordem

    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();

    cmesh->AutoBuild();


    return cmesh;
}
