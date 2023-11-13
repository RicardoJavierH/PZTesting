#include <TPZYSMPMatrix.h>
#include <pzstepsolver.h>
#include <fstream>
#include <iostream>

int main(){
    
    TPZFMatrix<STATE> matK;
    TPZFMatrix<STATE> rhs;
    matK.AutoFill(3, 3, SymProp::Sym);
    rhs.AutoFill(3, 1, SymProp::NonSym);

    //std::stringstream matname, rhsname;

    std::ofstream file("matK.txt");
    matK.Print("matK=",file,EMathematicaInput);
    std::ofstream file2("rhs.txt");
    rhs.Print("rhs=",file2,EMathematicaInput);
    file.close();
    
    //matK.Print(std::cout);

    TPZStepSolver<STATE> *solver = new TPZStepSolver<STATE>;
    solver->SetMatrix(&matK);
    solver->SetDirect(ELU);
    TPZFMatrix<STATE> solution;
    TPZFMatrix<STATE> residual;
    solver->Solve(rhs, solution, &residual);
    residual.Print(std::cout<<"residual:\n");
    
    solution.Print(std::cout << "solution:\n");
    std::ofstream outfile("solution.txt");
    solution.Print("sol=",outfile,EMathematicaInput);
    outfile.close();
    
    solution.Resize(3,1);
    TPZFMatrix<STATE> auxvec(3,1);
    matK.SetIsDecomposed(ENoDecompose);
    matK.Multiply(solution, auxvec);
    auxvec.Print(std::cout << "auxvec:\n");
    TPZFMatrix<STATE> result;
    auxvec.Subtract(rhs, result);
    result.Print(std::cout<< "result:\n");
    
    TPZFMatrix<STATE> residual2;
    matK.Residual(solution, rhs, residual2);
    residual2.Print(std::cout << "other_residual:\n");

    
//    if(0){
//        TPZFYsmpMatrix<STATE> sparse_mat;
//        sparse_mat.AutoFill(3, 3, SymProp::NonSym);
//
//        {
//            TPZMatrix<STATE> *mat_ptr_1 = &sparse_mat;
//            TPZMatrix<STATE> *mat_ptr_2 = sparse_mat.Clone();
//            mat_ptr_2->Zero();
//            *mat_ptr_2 = *mat_ptr_1;
//            for(int i = 0; i < 3; i++){
//                for(int j = 0; j < 3; j++){
//                    std::cout<<"i "<<i<<" j "<<j
//                    <<" ptr1 "<<mat_ptr_1->GetVal(i,j)
//                    <<" ptr2 "<<mat_ptr_2->GetVal(i,j)<<'\n';
//                }
//            }
//            delete mat_ptr_2;
//        }
//    }
        
    return 0;
}
