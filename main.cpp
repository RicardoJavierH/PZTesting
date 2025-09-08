#include <TPZYSMPMatrix.h>
#include <pzstepsolver.h>
#include <fstream>
#include <iostream>
#include <pzgmesh.h> //for TPZGeoMesh
#include "TPZVTKGeoMesh.h"
#include "TPZAnalyticSolution.h"
#include "DarcyFlow/TPZHybridDarcyFlow.h"//I can't invoke only the TPZDarcyFlow class

TPZGeoMesh* CreateTriangLShapeMesh(int nel, TPZVec<int>& bcids);
void UniformRefinement(int nDiv, TPZGeoMesh* gmesh);
TPZCompMesh* InsertCMeshH1(TPZGeoMesh* geomesh,TLaplaceExample1* exactsol, int intorder, int porder);


int main(){
    
    int porder = 1;
    int numinitialref = 0; // Number of refinements to be applied to the initial mesh
    int nthreads = 0;
    int integrationorder = 11;
    std::string topology = "Triangular"; //Triangular, Quadrilateral

    std::string problemname = "ESinMark";//ESinSin,ESinMark,EConst,EBubble2D,ESteepWave;
    TLaplaceExample1 aux, exact;
    exact.fExact = aux.ESinMark2;//ESinMark//ESinSin//ESinSinDirNonHom

    // Create geometric mesh
    TPZManVector<int, 8> Lshape_bcids(8, -1);
    TPZGeoMesh *gmesh = nullptr;
    int nelems= 6;
    gmesh = CreateTriangLShapeMesh(nelems, Lshape_bcids);
    //gmesh->Print();
    
    std::ofstream salida("mallageometrica.txt");
    gmesh->Print(salida);
    
    std::ofstream out("mallageom.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    
    UniformRefinement(numinitialref, gmesh);
    std::ofstream out2("mallageomrefinada.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2);
    
    //Create computational mesh
    TPZCompMesh* cmesh = nullptr;
    cmesh = InsertCMeshH1(gmesh, &exact, integrationorder, porder);
    std::ofstream out3("mallacomputac.txt");
    cmesh->Print(out3);
    return 0;
}

TPZGeoMesh* CreateTriangLShapeMesh(int nel, TPZVec<int>& bcids){
    
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
    
    // Creates quadrilateral element.
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


TPZCompMesh* InsertCMeshH1(TPZGeoMesh* geomesh,TLaplaceExample1* exactsol, int intorder, int porder) {

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
