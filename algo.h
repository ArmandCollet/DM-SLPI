#include <Dense>
#include<tuple>
#include<vector>
#include<iostream>
#include<Sparse>


//------------------------------ Gradient à Pas Optimal---------------------------------------

void GradPasOptimal(const Eigen::SparseMatrix<double> A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon,const int kmax, Eigen:: VectorXd & x);
void GradPasOptimal (const Eigen::SparseMatrix<double> A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x, std::string file);

//------------------------------Résidus minimum------------------------------------------------

void ResMin(const Eigen::SparseMatrix<double> A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x);
void ResMin (const Eigen::SparseMatrix<double> A,const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x, std::string file);

void ResMin_cond_gauche(Eigen::SparseMatrix<double> A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x);

void ResMin_cond_droite(Eigen::SparseMatrix<double> A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x);

void ResMin_cond_droite_flex(Eigen::SparseMatrix<double> A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x);


//-------------------------------GMRes--------------------------------------------------------

void GMRes(const Eigen::SparseMatrix<double> A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x, const int m);
void GMRes(const Eigen::SparseMatrix<double> A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x, const int m, std::string file);

void GMRes_preconditionne(const Eigen::SparseMatrix<double> A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x, const int m);

void GMRes_flex (const Eigen::SparseMatrix<double> A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x, const int m);


//--------------------------------------- Lecture des matrices---------------------------------

Eigen::SparseMatrix<double> Lecture_Matrice_A(std::string fichier);

Eigen::SparseMatrix<double> Lecture_Matrice_A_2(std::string fichier);

Eigen::VectorXd Lecture_Matrice_b(std::string fichier);





/*Eigen::MatrixXd ArnoldiH(const Eigen::MatrixXd A, Eigen::VectorXd & v);
  std::vector<Eigen::VectorXd> ArnoldiV(const Eigen::MatrixXd A, Eigen::VectorXd & v);*/

std::vector<Eigen::MatrixXd> Arnoldi(const Eigen::SparseMatrix<double> A, Eigen::VectorXd & v, const int m);

std::vector<Eigen::MatrixXd> Arnoldi_2(const Eigen::SparseMatrix<double> A,const Eigen::SparseMatrix<double> L,const Eigen::SparseMatrix<double> U, Eigen::VectorXd & v, const int m);

std::vector<Eigen::MatrixXd> Arnoldi_3(const Eigen::SparseMatrix<double> A, Eigen::VectorXd & v, const int m);

void GivensOpt(const Eigen::MatrixXd A, Eigen::MatrixXd & Q, Eigen::MatrixXd & R);

void resol_syst_triang_sup(const Eigen::MatrixXd A, Eigen::VectorXd & y, const Eigen::VectorXd b);

Eigen::VectorXd Resol_LU(Eigen::SparseMatrix<double> L, Eigen::SparseMatrix<double> U,Eigen::VectorXd b);
