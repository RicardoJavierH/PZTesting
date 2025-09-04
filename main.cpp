#include <TPZYSMPMatrix.h>
#include <pzstepsolver.h>
#include <fstream>
#include <iostream>
#include <pzgmesh.h> //for TPZGeoMesh
#include "TPZVTKGeoMesh.h"


TPZGeoMesh* CreateTriangLShapeMesh(int nel, TPZVec<int>& bcids);


int main(){
    TPZManVector<int, 8> Lshape_bcids(8, -1);
    TPZGeoMesh *gmesh;
    int nelems= 6;
    
    gmesh = CreateTriangLShapeMesh(nelems, Lshape_bcids);
    //gmesh->Print();
    
    std::ofstream salida("mallageometrica.txt");
    gmesh->Print(salida);
    
    std::ofstream out("mallarefinada.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
        
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
