#include <iostream>
#include <fstream>
#include "algo.h"
#include <vector>
#include <tuple>
using namespace std;
using namespace Eigen;

// -------------------------------------Gradient à Pas Optimal---------------------------------------------------

void GradPasOptimal (const SparseMatrix<double> A, const VectorXd b, const VectorXd x0, const double epsilon, const int kmax, VectorXd & x)
{
  //Initialisation
  int k=0;
  VectorXd r(b.size()), z(b.size());
  double alpha;
  r=b-A*x0;
  x=x0;

  //Boucle
  while((r.norm()>epsilon) && (k<=kmax))
  {
    z = A*r;
    alpha = r.dot(r) / z.dot(r);
    x = x + alpha*r;
    r = r - alpha*z;
    k+=1;
  }

  cout << "r="<<r.norm() << endl;
  cout<<"Nombre d'itérations ="<<k<<endl;

  if (k>kmax)
  {
    cout<<"Tolérance non atteinte: "<<endl;
  }

}

//sauvegarde les normes dans un fichier si on donne le nom du fichier
void GradPasOptimal (const SparseMatrix<double> A, const VectorXd b, const VectorXd x0, const double epsilon, const int kmax, VectorXd & x, std::string file)
{
  //Initialisation
  int k=0;
  VectorXd r(b.size()), z(b.size());
  double alpha;
  ofstream flux;
  flux.open(file);
  flux << "# itération k  ;  norme r"<<endl;
  r=b-A*x0;
  x=x0;
  //Boucle
  while((r.norm()>epsilon) && (k<=kmax))
  {
    z = A*r;
    alpha = r.dot(r) / z.dot(r);
    x = x + alpha*r;
    r = r - alpha*z;
    k+=1;
    flux << k << " " << r.norm() <<endl;
  }
  flux.close();
  cout << "r="<<r.norm() << endl;
  cout<<"Nombre d'itérations ="<<k<<endl;

  if (k>kmax)
  {
    cout<<"Tolérance non atteinte: "<<endl;
  }


}

//Algo identique à celui du prof
//(Test matrice aléatoire Q1)
//Les valeurs de x changent bien à chaque itération
//L'algo ne converge pas toujours en seulement 2 itérations.
//Testé avec des valeurs de x0 très différentes les unes des autres et nbr d'itérations varie aussi (atteint 9 avec ce que j'ai mis)
//Pour moi semble ok !


//-------------------------------------------------------Résidu minimium--------------------------------------------

void ResMin (const SparseMatrix<double> A,const VectorXd b, const VectorXd x0,const double epsilon,const int kmax, VectorXd & x)
{
  //Initialisation
  int k=0;
  VectorXd r(b.size()),z(b.size());
  double alpha;
  r=b-A*x0;
  x=x0;

  //Boucle
  while ((r.norm()>epsilon) && (k<=kmax))
  {
    z = A*r;
    alpha = r.dot(z) / z.dot(z);
    x = x + alpha*r;
    r = r - alpha*z;
    k+=1;
  }
  cout<<"r ="<<r.norm()<<endl;
  cout<<"Nombre d'itérations ="<<k<<endl;

  if (k>kmax)
  {
    cout<<"Tolérance non atteinte: "<<endl;
  }
}
//sauvegarde les normes dans un fichier si on donne le nom du fichier
void ResMin (const SparseMatrix<double> A,const VectorXd b, const VectorXd x0,const double epsilon,const int kmax, VectorXd & x, std::string file)
{
  //Initialisation
  int k=0;
  VectorXd r(b.size()),z(b.size());
  double alpha;
  ofstream flux;
  flux.open(file);
  flux << "# itération k  ;  norme r"<<endl;
  r=b-A*x0;
  x=x0;

  //Boucle
  while ((r.norm()>epsilon) && (k<=kmax))
  {
    z = A*r;
    alpha = r.dot(z) / z.dot(z);
    x = x + alpha*r;
    r = r - alpha*z;
    k+=1;
    flux << k << " " << r.norm() <<endl;
  }
  flux.close();
  cout<<"r ="<<r.norm()<<endl;
  cout<<"Nombre d'itérations ="<<k<<endl;

  if (k>kmax)
  {
    cout<<"Tolérance non atteinte: "<<endl;
  }
}

//Algo identique au prof
//(Test matrice aléatoire Q1)
//x varie bien à chaque itération
//Nombre d'itérations varie aussi en fonction du x0. Avec des valeurs wtf comme 175*i*i on monte à 8 itérations
//Algo me semble ok

//---------------------------------------Résolution LU-------------------------------------
//à voir si l'on s'en sert ou si l'on utilise la résolution traingulaire


//Résolution à partir d'une Décomposition LU
VectorXd Resol_LU(SparseMatrix <double> L, SparseMatrix<double> U, VectorXd b)
{
  int n=L.rows();
  VectorXd x(n), y(n);
  double s1, s2;

  y(0) = b(0);
  for(int i=1; i<n; i++)
  //for (int i=1; i<M.outerSize(); i++)
  {
    s1=0.;
    for(int k=0; k<i; k++)
    //for (SparseMatrix<double>::InnerIterator it(M,i); it; ++it)

    {
      //if(it.row()>it.col())
        //s1 += it.value()*y(it.row());
      s1 += L.coeffRef(i,k)*y(k);
    }
    y(i)=b(i)-s1;
  }
  x(n-1)=y(n-1)/U.coeffRef(n-1,n-1);


  for (int i=n-2; i>=0; i--)
  {
    s2=0.;

    for (int k=i+1;k<n;k++)
    //for (SparseMatrix<double>::InnerIterator it(M,i);it;++it)
    {
      //if (it.row()<it.col())
        //s2 += it.value()*x(it.row());
      s2 += U.coeffRef(i,k)*x(k);
    }
    x(i) = (y(i)-s2)/U.coeffRef(i,i);
  }
  return x;
}


//-----------------------------------------------Résidu minimum préconditionné à gauche--------------------------------------

void ResMin_cond_gauche(Eigen::SparseMatrix<double> A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x)
{

  // Initialisation //

  int m=A.rows();
  VectorXd r(b.size()),z(b.size()), z1(b.size()), q(m), q1(m), w(m);
  SparseMatrix<double> L(m,m), U(m,m), E(m,m), F(m,m), D(m,m), D_1(m,m);
  SparseLU<SparseMatrix<double> , COLAMDOrdering<int> > solver;
  double alpha;
  r=b-A*x0;

  // Création de M sous forme  de la décomposition LU //

  for (int i=0; i<A.outerSize(); ++i)
  {
    for (SparseMatrix<double>::InnerIterator it(A,i); it; ++it)
    {

      if (it.row()==it.col())
      {
        D.coeffRef(it.row(),it.col()) = it.value();
        D_1.coeffRef(it.row(),it.col()) = 1./it.value();
      }
      if (it.row()>it.col())
      {
        E.coeffRef(it.row(),it.col()) = -it.value();
      }

      if (it.row()<it.col())
      {
        F.coeffRef(it.row(),it.col()) = -it.value();
      }
    }
  }
  //cout<<"A="<<endl<<A<<endl;cout<<"D="<<endl<<D<<endl;cout<<"E="<<endl<<E<<endl;cout<<"F="<<endl<<F<<endl;
  L = (D-E)*D_1; U = (D-F);
  //cout<<L*U-(D-E)*D_1*(D-F)<<endl;
  x = x0;
  q1 = L.triangularView<Lower>().solve(r);
  q = U.triangularView<Upper>().solve(q1);
  //q = Resol_LU(L,U,r);

  int k=0;

  while ((r.norm()>epsilon) && (k<=kmax))
  {
    //cout<<"RM préconditionné à gauche, itération n°"<<k<<endl;
    w = A*q;
    //z = Resol_LU(L,U,w);
    z1 = L.triangularView<Lower>().solve(w);
    z = U.triangularView<Upper>().solve(z1);
    alpha = q.dot(z)/z.dot(z);
    x = x + alpha*q;
    r = r - alpha*w;
    q = q - alpha*z;
    k += 1;
  }
  cout<<"r ="<<r.norm()<<endl;
  cout<<"Nombre d'itérations ="<<k<<endl;
  if (k>kmax)
  {
    cout<<"Tolérance non atteinte: "<<endl;
  }
}

//sauvegarde de la solution
void ResMin_cond_gauche(Eigen::SparseMatrix<double> A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x, std::string file)
{

  // Initialisation //

  int m=A.rows();
  VectorXd r(b.size()),z(b.size()), z1(b.size()), q(m), q1(m), w(m);
  SparseMatrix<double> L(m,m), U(m,m), E(m,m), F(m,m), D(m,m), D_1(m,m);
  SparseLU<SparseMatrix<double> , COLAMDOrdering<int> > solver;
  double alpha;
  ofstream flux;
  flux.open(file);
  flux << "# itération k  ;  norme r"<<endl;
  r=b-A*x0;

  // Création de M sous forme  de la décomposition LU //

  for (int i=0; i<A.outerSize(); ++i)
  {
    for (SparseMatrix<double>::InnerIterator it(A,i); it; ++it)
    {

      if (it.row()==it.col())
      {
        D.coeffRef(it.row(),it.col()) = it.value();
        D_1.coeffRef(it.row(),it.col()) = 1./it.value();
      }
      if (it.row()>it.col())
      {
        E.coeffRef(it.row(),it.col()) = -it.value();
      }

      if (it.row()<it.col())
      {
        F.coeffRef(it.row(),it.col()) = -it.value();
      }
    }
  }
  //cout<<"A="<<endl<<A<<endl;cout<<"D="<<endl<<D<<endl;cout<<"E="<<endl<<E<<endl;cout<<"F="<<endl<<F<<endl;
  L = (D-E)*D_1; U = (D-F);
  //cout<<L*U-(D-E)*D_1*(D-F)<<endl;
  x = x0;
  q1 = L.triangularView<Lower>().solve(r);
  q = U.triangularView<Upper>().solve(q1);
  //q = Resol_LU(L,U,r);

  int k=0;

  while ((r.norm()>epsilon) && (k<=kmax))
  {
    //cout<<"RM préconditionné à gauche, itération n°"<<k<<endl;
    w = A*q;
    //z = Resol_LU(L,U,w);
    z1 = L.triangularView<Lower>().solve(w);
    z = U.triangularView<Upper>().solve(z1);
    alpha = q.dot(z)/z.dot(z);
    x = x + alpha*q;
    r = r - alpha*w;
    q = q - alpha*z;
    k += 1;
    flux << k << " " << r.norm() <<endl;
  }
  flux.close();

  cout<<"r ="<<r.norm()<<endl;
  cout<<"Nombre d'itérations ="<<k<<endl;
  if (k>kmax)
  {
    cout<<"Tolérance non atteinte: "<<endl;
  }
}


//------------------------------------------------Résidu minimum conditionné à droite-------------------------------------


void ResMin_cond_droite(Eigen::SparseMatrix<double> A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x)
{
  int m=A.rows();
  VectorXd r(b.size()); r=b-A*x0;
  VectorXd z(b.size()), z1(b.size());
  VectorXd q(m);
  VectorXd w(m);
  double alpha;
  SparseMatrix<double> L(m,m), U(m,m), E(m,m), F(m,m), D(m,m), D_1(m,m);

  for (int i=0; i<A.outerSize(); ++i)
  {
    for (SparseMatrix<double>::InnerIterator it(A,i); it; ++it)
    {

      if (it.row()==it.col())
      {
        D.coeffRef(it.row(),it.col()) = it.value();
        D_1.coeffRef(it.row(),it.col()) = 1./it.value();
      }
      if (it.row()>it.col())
      {
        E.coeffRef(it.row(),it.col()) = -it.value();
      }

      if (it.row()<it.col())
      {
        F.coeffRef(it.row(),it.col()) = -it.value();
      }
    }
  }

  L = (D-E)*D_1; U = (D-F);

  x = x0;

  int k=0;

  while ((r.norm()>epsilon) && (k<=kmax))
  {
    //cout<<"RM préconditionné à gauche, itération n°"<<k<<endl;

    //z = Resol_LU(L,U,r);
    z1 = L.triangularView<Lower>().solve(r);
    z = U.triangularView<Upper>().solve(z1);
    w = A*z;
    alpha = r.dot(w)/w.dot(w);
    x = x + alpha*z;
    r = r - alpha*w;

    k += 1;
  }
  cout<<"r ="<<r.norm()<<endl;
  cout<<"Nombre d'itérations ="<<k<<endl;
  if (k>kmax)
  {
    cout<<"Tolérance non atteinte: "<<endl;
  }
}

//sauvegarde de solution
void ResMin_cond_droite(Eigen::SparseMatrix<double> A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x, std::string file)
{
  int m=A.rows();
  VectorXd r(b.size()); r=b-A*x0;
  VectorXd z(b.size()), z1(b.size());
  VectorXd q(m);
  VectorXd w(m);
  double alpha;
  SparseMatrix<double> L(m,m), U(m,m), E(m,m), F(m,m), D(m,m), D_1(m,m);
   ofstream flux;
  flux.open(file);
  flux << "# itération k  ;  norme r"<<endl;
  for (int i=0; i<A.outerSize(); ++i)
  {
    for (SparseMatrix<double>::InnerIterator it(A,i); it; ++it)
    {

      if (it.row()==it.col())
      {
        D.coeffRef(it.row(),it.col()) = it.value();
        D_1.coeffRef(it.row(),it.col()) = 1./it.value();
      }
      if (it.row()>it.col())
      {
        E.coeffRef(it.row(),it.col()) = -it.value();
      }

      if (it.row()<it.col())
      {
        F.coeffRef(it.row(),it.col()) = -it.value();
      }
    }
  }

  L = (D-E)*D_1; U = (D-F);

  x = x0;

  int k=0;

  while ((r.norm()>epsilon) && (k<=kmax))
  {
    //cout<<"RM préconditionné à gauche, itération n°"<<k<<endl;

    //z = Resol_LU(L,U,r);
    z1 = L.triangularView<Lower>().solve(r);
    z = U.triangularView<Upper>().solve(z1);
    w = A*z;
    alpha = r.dot(w)/w.dot(w);
    x = x + alpha*z;
    r = r - alpha*w;

    k += 1;
    flux << k << " " << r.norm() <<endl;
  }
  flux.close();

  cout<<"r ="<<r.norm()<<endl;
  cout<<"Nombre d'itérations ="<<k<<endl;
  if (k>kmax)
  {
    cout<<"Tolérance non atteinte: "<<endl;
  }
}

//------------------------Résidu minimum préconditionné à droite flex----------------------------

void ResMin_cond_droite_flex(Eigen::SparseMatrix<double> A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x)
{
  int m=A.rows();
  VectorXd r(b.size()); r=b-A*x0;
  VectorXd z(b.size()), z1(b.size()), z0(b.size());
  VectorXd q(m);
  VectorXd w(m);
  double alpha;
  SparseMatrix<double> L(m,m), U(m,m), E(m,m), F(m,m), D(m,m), D_1(m,m);

  x = x0;

  int k=0;
  z0.setZero();
  while ((r.norm()>epsilon) && (k<=kmax))
  {
    //cout<<"RM préconditionné à gauche, itération n°"<<k<<endl;

    //z = Resol_LU(L,U,r);
    //z1 = L.triangularView<Lower>().solve(r);
    //z = U.triangularView<Upper>().solve(z1);
    ResMin (A,r,z0,epsilon,10,z);
    z0=z;
    w = A*z;
    alpha = r.dot(w)/w.dot(w);
    x = x + alpha*z;
    r = r - alpha*w;

    k += 1;
  }
  cout<<"r ="<<r.norm()<<endl;
  cout<<"Nombre d'itérations ="<<k<<endl;
  if (k>kmax)
  {
    cout<<"Tolérance non atteinte: "<<endl;
  }
}

//sauvegarde solution
void ResMin_cond_droite_flex(Eigen::SparseMatrix<double> A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x, std::string file)
{
  int m=A.rows();
  VectorXd r(b.size()); r=b-A*x0;
  VectorXd z(b.size()), z1(b.size()), z0(b.size());
  VectorXd q(m);
  VectorXd w(m);
  double alpha;
  SparseMatrix<double> L(m,m), U(m,m), E(m,m), F(m,m), D(m,m), D_1(m,m);


   ofstream flux;
  flux.open(file);
  flux << "# itération k  ;  norme r"<<endl;
  x = x0;

  int k=0;
  z0.setZero();
  while ((r.norm()>epsilon) && (k<=kmax))
  {
    //cout<<"RM préconditionné à gauche, itération n°"<<k<<endl;

    //z = Resol_LU(L,U,r);
    //z1 = L.triangularView<Lower>().solve(r);
    //z = U.triangularView<Upper>().solve(z1);
    ResMin (A,r,z0,epsilon,10,z);
    z0=z;
    w = A*z;
    alpha = r.dot(w)/w.dot(w);
    x = x + alpha*z;
    r = r - alpha*w;

    k += 1;
    flux << k << " " << r.norm() <<endl;
  }
  flux.close();

  cout<<"r ="<<r.norm()<<endl;
  cout<<"Nombre d'itérations ="<<k<<endl;
  if (k>kmax)
  {
    cout<<"Tolérance non atteinte: "<<endl;
  }
}

//-------------------------------------GMRes-----------------------------------------------------



//Arnoldi
std::vector<Eigen::MatrixXd> Arnoldi(const Eigen::SparseMatrix<double> A, Eigen::VectorXd & v, const int m)
{
  //Déclaration des variables

  MatrixXd Vm_1(v.size(),m+1); Vm_1 = MatrixXd::Constant(v.size(),m+1,0.);
  MatrixXd Vm(v.size(),m);
  MatrixXd Hm(m,m);
  MatrixXd Hm_barre(m+1,m); Hm_barre = MatrixXd::Constant(m+1,m,0.);
  MatrixXd w(v.size(),m);
  vector<Eigen::MatrixXd> HmVm;


  //Initialisation
  Vm_1.col(0)=v/v.norm();

  //Boucles
  for (int j=0; j<m;j++)
  {
    w.col(j)=A*Vm_1.col(j);
    for (int i=0;i<=j;i++)
    {
      Hm_barre(i,j)=w.col(j).dot(Vm_1.col(i));
      w.col(j)=w.col(j)-Hm_barre(i,j)*Vm_1.col(i);
    }
    Hm_barre(j+1,j)=(w.col(j)).norm();
    if (Hm_barre(j+1,j)==0.)
    break;
    Vm_1.col(j+1)=w.col(j)/Hm_barre(j+1,j);

  }

  //Extraction de Vm et Hm à partie de Vm+1 et Hm_barre
  for(int i=0;i<m;i++)
  {
    Vm.col(i)=Vm_1.col(i);Hm.row(i)=Hm_barre.row(i);
  }

  //Remplissage du vecteur contenant Hm, Hm_barre, Vm et Vm+1
  HmVm.push_back(Hm);
  HmVm.push_back(Hm_barre);
  HmVm.push_back(Vm);
  HmVm.push_back(Vm_1);

  return HmVm;
}

//Givens honnête
void GivensOpt(const Eigen::MatrixXd A, Eigen::MatrixXd & Q, Eigen::MatrixXd & R)
{
  int m; m = A.rows(); int n = A.cols();
  R = A ; Q = MatrixXd::Identity(m,m);

  for (int i=0;i<m;i++)
  {

    for (int j=0; j<i;j++)
    {
      double a = sqrt( pow(R(j,j),2) + pow(R(i,j),2) );// wj_intermediaire = L.triangularView<Lower>().solve(Vm_1.col(j));
    // wj0 = U.triangularView<Upper>().solve(wj_intermediaire);
      double c = R(j,j)/a;
      double s = -R(i,j)/a;

      for (int l=0; l<m ; l++)
      {
        if (l<n)
        {
          double Ril_inter = c*R(i,l) + s*R(j,l);
          R(j,l) = c*R(j,l) - s*R(i,l);
          R(i,l) = Ril_inter;
        }
      }

      for (int l=0;l<m;l++)
      {
          double Qil_inter = s*Q(l,j)+c*Q(l,i);
          Q(l,j) = -s*Q(l,i)+c*Q(l,j);
          Q(l,i) = Qil_inter;
      }
    }
  }
}

//Résolution système triangulaire supérieur
void resol_syst_triang_sup(const Eigen::MatrixXd A, Eigen::VectorXd & y, const Eigen::VectorXd b)
{
  double S;
  int m; m = A.rows();

  for (int i=m-1; i>=0;i--)
  {
    S = 0.;
    for (int j=m-1;j>i;j--)
    {
      S += A(i,j)*y(j);
    }
    y(i) = (b(i)-S)/A(i,i);
  }

}

//GMRes
void GMRes(const Eigen::SparseMatrix<double> A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x, const int m)
{
  //Initialisation
  int k=0;
  VectorXd r(b.size()); r=b-A*x0;
  VectorXd y(m);
  double beta; beta=r.norm(); double gamma;
  MatrixXd Hm(m,m); MatrixXd Vm(b.size(),m);MatrixXd Hm_barre(m+1,m);MatrixXd Vm_1(b.size(),m+1);
  vector<Eigen::MatrixXd> HmVm;
  Eigen::MatrixXd Qm(m,m); Eigen::MatrixXd Rm(m,m); Eigen::MatrixXd Qm_barre(m+1,m+1);Eigen::MatrixXd Rm_barre(m+1,m);
  Eigen::VectorXd gm_barre(m+1);Eigen::VectorXd gm(m);

  x=x0;

  while ((beta>epsilon) && (k<=kmax))
  {
    HmVm = Arnoldi(A,r,m);
    Hm = HmVm[0]; Vm = HmVm[2]; Hm_barre = HmVm[1]; Vm_1 = HmVm[3];
    //cout<<"VmT*Vm="<<endl;
    //cout<<Vm.transpose()*Vm<<endl;

    //Givens(Hm,Qm,Rm);
    GivensOpt(Hm_barre,Qm_barre,Rm_barre);

    gm_barre = beta*(Qm_barre.transpose()).col(0);
    gamma = gm_barre[m];

    for(int i=0;i<m;i++)
    {
      Rm.row(i) = Rm_barre.row(i); gm.row(i) = gm_barre.row(i);
    }

    resol_syst_triang_sup(Rm, y, gm);

    x += Vm * y;

    r = b-A*x;
    //r = gamma*Vm_1*((Qm_barre.transpose()).col(m));

    //beta = abs(gamma);
    beta = r.norm();
    cout<<"r ="<<r.norm()<<endl;
    cout<<"Nombre d'itérations ="<<k<<endl;
    k+=1;
  }
  cout<<"r ="<<r.norm()<<endl;
  cout<<"Nombre d'itérations ="<<k<<endl;
  if (k>kmax)
  {
    cout<<"Tolérance non atteinte: "<<endl;
  }

}

//sauvegarde les normes dans un fichier si on donne le nom du fichier
void GMRes(const Eigen::SparseMatrix<double> A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x, const int m, std::string file)
{
  //Initialisation
  int k=0;
  VectorXd r(b.size()); r=b-A*x0;
  VectorXd y(m);
  double beta; beta=r.norm(); double gamma;
  MatrixXd Hm(m,m); MatrixXd Vm(b.size(),m);MatrixXd Hm_barre(m+1,m);MatrixXd Vm_1(b.size(),m+1);
  vector<Eigen::MatrixXd> HmVm;
  Eigen::MatrixXd Qm(m,m); Eigen::MatrixXd Rm(m,m); Eigen::MatrixXd Qm_barre(m+1,m+1);Eigen::MatrixXd Rm_barre(m+1,m);
  Eigen::VectorXd gm_barre(m+1);Eigen::VectorXd gm(m);

  ofstream flux;
  flux.open(file);
  flux << "# itération k  ;  norme r"<<endl;

  x=x0;

  while ((beta>epsilon) && (k<=kmax))
  {
    HmVm = Arnoldi(A,r,m);
    Hm = HmVm[0]; Vm = HmVm[2]; Hm_barre = HmVm[1]; Vm_1 = HmVm[3];
    //cout<<"VmT*Vm="<<endl;
    //cout<<Vm.transpose()*Vm<<endl;

    //Givens(Hm,Qm,Rm);
    GivensOpt(Hm_barre,Qm_barre,Rm_barre);

    gm_barre = beta*(Qm_barre.transpose()).col(0);
    gamma = gm_barre[m];

    for(int i=0;i<m;i++)
    {
      Rm.row(i) = Rm_barre.row(i); gm.row(i) = gm_barre.row(i);
    }

    resol_syst_triang_sup(Rm, y, gm);

    x += Vm * y;

    r = b-A*x;
    //r = gamma*Vm_1*((Qm_barre.transpose()).col(m));

    //beta = abs(gamma);
    beta = r.norm();
    k+=1;
    flux << k << " " << r.norm() <<endl;
  }
  flux.close();
  cout<<"r ="<<r.norm()<<endl;
  cout<<"Nombre d'itérations ="<<k<<endl;
  if (k>kmax)
  {
    cout<<"Tolérance non atteinte: "<<endl;
  }

}

//---------------------------------GMRes préconditionné-------------------------------------------------//
std::vector<Eigen::MatrixXd> Arnoldi_2(const Eigen::SparseMatrix<double> A,const Eigen::SparseMatrix<double> L,const Eigen::SparseMatrix<double> U, Eigen::VectorXd & v, const int m)
{
  //Déclaration des variables

  MatrixXd Vm_1(v.size(),m+1); Vm_1 = MatrixXd::Constant(v.size(),m+1,0.);
  MatrixXd Vm(v.size(),m);   MatrixXd Zm(v.size(),m);
  MatrixXd Hm(m,m);
  MatrixXd Hm_barre(m+1,m); Hm_barre = MatrixXd::Constant(m+1,m,0.);
  MatrixXd w(v.size(),m);VectorXd wj0(w.rows()), wj_intermediaire(w.rows());
  vector<Eigen::MatrixXd> HmVm;




  //Initialisation
  Vm_1.col(0)=v/v.norm();

  //Boucles
  for (int j=0; j<m;j++)
  {
    wj_intermediaire = L.triangularView<Lower>().solve(Vm_1.col(j));
    wj0 = U.triangularView<Upper>().solve(wj_intermediaire);
    w.col(j)=A*wj0;
    Zm.col(j)=wj0;
    //w.col(j)=A*Vm_1.col(j);
    for (int i=0;i<=j;i++)
    {
      Hm_barre(i,j)=w.col(j).dot(Vm_1.col(i));
      w.col(j)=w.col(j)-Hm_barre(i,j)*Vm_1.col(i);
    }
    Hm_barre(j+1,j)=(w.col(j)).norm();
    if (Hm_barre(j+1,j)==0.)
    break;
    Vm_1.col(j+1)=w.col(j)/Hm_barre(j+1,j);

  }

  //Extraction de Vm et Hm à partie de Vm+1 et Hm_barre
  for(int i=0;i<m;i++)
  {
    Vm.col(i)=Vm_1.col(i);Hm.row(i)=Hm_barre.row(i);
  }
  //cout<<Vm.transpose()*Vm<<endl;
  //Remplissage du vecteur contenant Hm, Hm_barre, Vm et Vm+1
  HmVm.push_back(Hm);
  HmVm.push_back(Hm_barre);
  HmVm.push_back(Vm);
  HmVm.push_back(Vm_1);
  HmVm.push_back(Zm);

  return HmVm;
}

void GMRes_preconditionne(const Eigen::SparseMatrix<double> A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x, const int m)
{
  //Initialisation
  int k=0;
  VectorXd r(b.size()); r=b-A*x0;
  VectorXd y(m); VectorXd u(x0.size()), x_intermediaire(x0.size());
  double beta; beta=r.norm(); double gamma;
  MatrixXd Hm(m,m); MatrixXd Vm(b.size(),m);MatrixXd Hm_barre(m+1,m);MatrixXd Vm_1(b.size(),m+1); MatrixXd Zm(b.size(),m);
  vector<Eigen::MatrixXd> HmVm;
  Eigen::MatrixXd Qm(m,m); Eigen::MatrixXd Rm(m,m); Eigen::MatrixXd Qm_barre(m+1,m+1);Eigen::MatrixXd Rm_barre(m+1,m);
  Eigen::VectorXd gm_barre(m+1);Eigen::VectorXd gm(m);
  SparseMatrix<double> L(A.rows(),A.rows()), U(A.rows(),A.rows()), E(A.rows(),A.rows()), F(A.rows(),A.rows()), D(A.rows(),A.rows()), D_1(A.rows(),A.rows());

  for (int i=0; i<A.outerSize(); ++i)
  {
    for (SparseMatrix<double>::InnerIterator it(A,i); it; ++it)
    {

      if (it.row()==it.col())
      {
        D.coeffRef(it.row(),it.col()) = it.value();
        D_1.coeffRef(it.row(),it.col()) = 1./it.value();
      }
      if (it.row()>it.col())
      {
        E.coeffRef(it.row(),it.col()) = -it.value();
      }

      if (it.row()<it.col())
      {
        F.coeffRef(it.row(),it.col()) = -it.value();
      }
    }
  }

  L = (D-E)*D_1; U = (D-F);

  x=x0;

  while ((beta>epsilon) && (k<=kmax))
  {
    HmVm = Arnoldi_2(A,L,U,r,m);
    Hm = HmVm[0]; Vm = HmVm[2]; Hm_barre = HmVm[1]; Vm_1 = HmVm[3]; Zm = HmVm[4];
    //cout<<"VmT*Vm="<<endl;
    //cout<<Vm.transpose()*Vm<<endl;

    //Givens(Hm,Qm,Rm);
    GivensOpt(Hm_barre,Qm_barre,Rm_barre);

    gm_barre = beta*(Qm_barre.transpose()).col(0);
    //gamma = gm_barre[m];

    for(int i=0;i<m;i++)
    {
      Rm.row(i) = Rm_barre.row(i); gm.row(i) = gm_barre.row(i);
    }

    resol_syst_triang_sup(Rm, y, gm);

    x += Zm*y;

    r = b-A*x;
    //r = gamma*Vm_1*((Qm_barre.transpose()).col(m));

    //beta = abs(gamma);
    beta = r.norm();

    k+=1;
  }

  cout<<"r ="<<beta<<endl;
  cout<<"Nombre d'itérations ="<<k<<endl;
  if (k>kmax)
  {
    cout<<"Tolérance non atteinte: "<<endl;
  }

}
//sauvegarde solution
void GMRes_preconditionne(const Eigen::SparseMatrix<double> A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x, const int m, std::string file)
{
  //Initialisation
  int k=0;
  VectorXd r(b.size()); r=b-A*x0;
  VectorXd y(m); VectorXd u(x0.size()), x_intermediaire(x0.size());
  double beta; beta=r.norm(); double gamma;
  MatrixXd Hm(m,m); MatrixXd Vm(b.size(),m);MatrixXd Hm_barre(m+1,m);MatrixXd Vm_1(b.size(),m+1); MatrixXd Zm(b.size(),m);
  vector<Eigen::MatrixXd> HmVm;
  Eigen::MatrixXd Qm(m,m); Eigen::MatrixXd Rm(m,m); Eigen::MatrixXd Qm_barre(m+1,m+1);Eigen::MatrixXd Rm_barre(m+1,m);
  Eigen::VectorXd gm_barre(m+1);Eigen::VectorXd gm(m);
  SparseMatrix<double> L(A.rows(),A.rows()), U(A.rows(),A.rows()), E(A.rows(),A.rows()), F(A.rows(),A.rows()), D(A.rows(),A.rows()), D_1(A.rows(),A.rows());

  ofstream flux;
  flux.open(file);
  flux << "# itération k  ;  norme r"<<endl;

  for (int i=0; i<A.outerSize(); ++i)
  {
    for (SparseMatrix<double>::InnerIterator it(A,i); it; ++it)
    {

      if (it.row()==it.col())
      {
        D.coeffRef(it.row(),it.col()) = it.value();
        D_1.coeffRef(it.row(),it.col()) = 1./it.value();
      }
      if (it.row()>it.col())
      {
        E.coeffRef(it.row(),it.col()) = -it.value();
      }

      if (it.row()<it.col())
      {
        F.coeffRef(it.row(),it.col()) = -it.value();
      }
    }
  }

  L = (D-E)*D_1; U = (D-F);

  x=x0;

  while ((beta>epsilon) && (k<=kmax))
  {
    HmVm = Arnoldi_2(A,L,U,r,m);
    Hm = HmVm[0]; Vm = HmVm[2]; Hm_barre = HmVm[1]; Vm_1 = HmVm[3]; Zm = HmVm[4];
    //cout<<"VmT*Vm="<<endl;
    //cout<<Vm.transpose()*Vm<<endl;

    //Givens(Hm,Qm,Rm);
    GivensOpt(Hm_barre,Qm_barre,Rm_barre);

    gm_barre = beta*(Qm_barre.transpose()).col(0);
    //gamma = gm_barre[m];

    for(int i=0;i<m;i++)
    {
      Rm.row(i) = Rm_barre.row(i); gm.row(i) = gm_barre.row(i);
    }

    resol_syst_triang_sup(Rm, y, gm);

    x += Zm*y;

    r = b-A*x;
    //r = gamma*Vm_1*((Qm_barre.transpose()).col(m));

    //beta = abs(gamma);
    beta = r.norm();

    k+=1;
    flux << k << " " << r.norm() <<endl;
  }
  flux.close();


  cout<<"r ="<<beta<<endl;
  cout<<"Nombre d'itérations ="<<k<<endl;
  if (k>kmax)
  {
    cout<<"Tolérance non atteinte: "<<endl;
  }

}
//----------------------------------------Flex GMRes standard--------------------------------------------------------//
std::vector<Eigen::MatrixXd> Arnoldi_3(const Eigen::SparseMatrix<double> A, Eigen::VectorXd & v, const int m)
{
  //Déclaration des variables

  MatrixXd Vm_1(v.size(),m+1); Vm_1 = MatrixXd::Constant(v.size(),m+1,0.);
  MatrixXd Vm(v.size(),m);   MatrixXd Zm(v.size(),m);
  MatrixXd Hm(m,m);
  MatrixXd Hm_barre(m+1,m); Hm_barre = MatrixXd::Constant(m+1,m,0.);
  MatrixXd w(v.size(),m);VectorXd wj0(w.rows()), wj_intermediaire(w.rows()), wj(w.rows());
  vector<Eigen::MatrixXd> HmVm;




  //Initialisation
  Vm_1.col(0)=v/v.norm();

  //Boucles
  for (int j=0; j<m;j++)
  {
    GMRes (A, Vm_1.col(j), wj0, 0.000001, 10, wj, 5);
    w.col(j)=A*wj;
    wj0=wj;
    Zm.col(j)=wj0;
    //w.col(j)=A*Vm_1.col(j);
    for (int i=0;i<=j;i++)
    {
      Hm_barre(i,j)=w.col(j).dot(Vm_1.col(i));
      w.col(j)=w.col(j)-Hm_barre(i,j)*Vm_1.col(i);
    }
    Hm_barre(j+1,j)=(w.col(j)).norm();
    if (Hm_barre(j+1,j)==0.)
    break;
    Vm_1.col(j+1)=w.col(j)/Hm_barre(j+1,j);

  }

  //Extraction de Vm et Hm à partie de Vm+1 et Hm_barre
  for(int i=0;i<m;i++)
  {
    Vm.col(i)=Vm_1.col(i);Hm.row(i)=Hm_barre.row(i);
  }
  //cout<<Vm.transpose()*Vm<<endl;
  //Remplissage du vecteur contenant Hm, Hm_barre, Vm et Vm+1
  HmVm.push_back(Hm);
  HmVm.push_back(Hm_barre);
  HmVm.push_back(Vm);
  HmVm.push_back(Vm_1);
  HmVm.push_back(Zm);

  return HmVm;
}

void GMRes_flex(const Eigen::SparseMatrix<double> A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x, const int m)
{
  //Initialisation
  int k=0;
  VectorXd r(b.size()); r=b-A*x0;
  VectorXd y(m); VectorXd u(x0.size()), x_intermediaire(x0.size());
  double beta; beta=r.norm(); double gamma;
  MatrixXd Hm(m,m); MatrixXd Vm(b.size(),m);MatrixXd Hm_barre(m+1,m);MatrixXd Vm_1(b.size(),m+1); MatrixXd Zm(b.size(),m);
  vector<Eigen::MatrixXd> HmVm;
  Eigen::MatrixXd Qm(m,m); Eigen::MatrixXd Rm(m,m); Eigen::MatrixXd Qm_barre(m+1,m+1);Eigen::MatrixXd Rm_barre(m+1,m);
  Eigen::VectorXd gm_barre(m+1);Eigen::VectorXd gm(m);
  SparseMatrix<double> L(A.rows(),A.rows()), U(A.rows(),A.rows()), E(A.rows(),A.rows()), F(A.rows(),A.rows()), D(A.rows(),A.rows()), D_1(A.rows(),A.rows());


  x=x0;

  while ((beta>epsilon) && (k<=kmax))
  {
    HmVm = Arnoldi_3(A,r,m);
    Hm = HmVm[0]; Vm = HmVm[2]; Hm_barre = HmVm[1]; Vm_1 = HmVm[3]; Zm = HmVm[4];
    //cout<<"VmT*Vm="<<endl;
    //cout<<Vm.transpose()*Vm<<endl;

    //Givens(Hm,Qm,Rm);
    GivensOpt(Hm_barre,Qm_barre,Rm_barre);

    gm_barre = beta*(Qm_barre.transpose()).col(0);
    //gamma = gm_barre[m];

    for(int i=0;i<m;i++)
    {
      Rm.row(i) = Rm_barre.row(i); gm.row(i) = gm_barre.row(i);
    }

    resol_syst_triang_sup(Rm, y, gm);

    x += Zm*y;

    r = b-A*x;

    beta = r.norm();cout<<"Iteration GMRes flex"<<beta<<endl;

    k+=1;
  }

  cout<<"r ="<<beta<<endl;
  cout<<"Nombre d'itérations ="<<k<<endl;
  if (k>kmax)
  {
    cout<<"Tolérance non atteinte: "<<endl;
  }

}

//save solution
void GMRes_flex(const Eigen::SparseMatrix<double> A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x, const int m, std::string file)
{
  //Initialisation
  int k=0;
  VectorXd r(b.size()); r=b-A*x0;
  VectorXd y(m); VectorXd u(x0.size()), x_intermediaire(x0.size());
  double beta; beta=r.norm(); double gamma;
  MatrixXd Hm(m,m); MatrixXd Vm(b.size(),m);MatrixXd Hm_barre(m+1,m);MatrixXd Vm_1(b.size(),m+1); MatrixXd Zm(b.size(),m);
  vector<Eigen::MatrixXd> HmVm;
  Eigen::MatrixXd Qm(m,m); Eigen::MatrixXd Rm(m,m); Eigen::MatrixXd Qm_barre(m+1,m+1);Eigen::MatrixXd Rm_barre(m+1,m);
  Eigen::VectorXd gm_barre(m+1);Eigen::VectorXd gm(m);
  SparseMatrix<double> L(A.rows(),A.rows()), U(A.rows(),A.rows()), E(A.rows(),A.rows()), F(A.rows(),A.rows()), D(A.rows(),A.rows()), D_1(A.rows(),A.rows());

  ofstream flux;
  flux.open(file);
  flux << "# itération k  ;  norme r"<<endl;
  x=x0;

  while ((beta>epsilon) && (k<=kmax))
  {
    HmVm = Arnoldi_3(A,r,m);
    Hm = HmVm[0]; Vm = HmVm[2]; Hm_barre = HmVm[1]; Vm_1 = HmVm[3]; Zm = HmVm[4];
    //cout<<"VmT*Vm="<<endl;
    //cout<<Vm.transpose()*Vm<<endl;

    //Givens(Hm,Qm,Rm);
    GivensOpt(Hm_barre,Qm_barre,Rm_barre);

    gm_barre = beta*(Qm_barre.transpose()).col(0);
    //gamma = gm_barre[m];

    for(int i=0;i<m;i++)
    {
      Rm.row(i) = Rm_barre.row(i); gm.row(i) = gm_barre.row(i);
    }

    resol_syst_triang_sup(Rm, y, gm);

    x += Zm*y;

    r = b-A*x;

    beta = r.norm();cout<<"Iteration GMRes flex"<<beta<<endl;

    k+=1;
    flux << k << " " << r.norm() <<endl;
  }
  flux.close();
  cout<<"r ="<<beta<<endl;
  cout<<"Nombre d'itérations ="<<k<<endl;
  if (k>kmax)
  {
    cout<<"Tolérance non atteinte: "<<endl;
  }

}
//------------------------------------------LECTURE MATRICES-------------------------------------------//

SparseMatrix <double>  Lecture_Matrice_A (string fichier)
{
  cout << "Lecture du fichier : "<<fichier<<endl;
  ifstream flux_A(fichier);
  SparseMatrix <double> A;
  int n,m,kmax;
  string ligne;
  if (flux_A)
    {
      for (int i=0; i<34 ; i++)
	{
	  getline(flux_A,ligne);
	}
      flux_A >> n;
      flux_A >> m;
      flux_A >> kmax;
      cout << "n, m, imax : "<<n<<" "<<m<<" "<<kmax<<endl<<endl;
      A.resize(n,m);
      vector<Triplet<double>> triplets;
      double a(0);
      int i(0);
      int j(0);
      for (int k=0; k<kmax; k++)
	{

	  flux_A >> i >> j >> a ;

    if (i!=j)
    {
	     triplets.push_back({i-1,j-1,a});
	     triplets.push_back({j-1,i-1,a});
     }

     else
     {
       triplets.push_back({i-1,j-1,a});
     }
	}
      A.setFromTriplets(triplets.begin(), triplets.end());
    }
  else
    {
      cout << "Ce fichier n'existe pas"<<endl;
    }
  flux_A.close();
  // cout << "A(1,1) = "<<A.coeffRef(0,0)<<endl
  // <<A.coeffRef(19,0)<<" "<<A.coeffRef(0,19)<<endl;
  // cout<< "dernier terme diagonale : "<<A.coeffRef(25709,25709)<<endl
  // <<A.coeffRef(25709,25708)<<endl;
  return A;
}

VectorXd Lecture_Matrice_b(string fichier)
{
  cout << "Lecture du fichier : "<<fichier<<endl;
  ifstream flux_b(fichier); //../smt/smt_b.mtx
  VectorXd b;
  int kmax;
  double a;
  string ligne;
  if (flux_b)
    {
      for (int i=0; i<6 ; i++)
	{
	  getline(flux_b,ligne);
	}
      flux_b >> kmax >> a; //a va être réécrit
      cout << "kmax : "<<kmax<<endl<<endl;
      b.resize(kmax);
      for (int k=1; k<kmax+1 ; k++)
	{
	  flux_b >> a;
	  b[k-1] = a;
	}
    }
  else
    {
      cout<<"Ce fichier n'existe pas"<<endl;
    }
  flux_b.close();
  return b;
}

SparseMatrix <double>  Lecture_Matrice_A_2 (string fichier)
{
  cout << "Lecture du fichier : "<<fichier<<endl;
  ifstream flux_A(fichier);
  SparseMatrix <double> A;
  int n,m,kmax;
  string ligne;
  if (flux_A)
    {
      for (int i=0; i<27 ; i++)
	{
	  getline(flux_A,ligne);
	}
      flux_A >> n;
      flux_A >> m;
      flux_A >> kmax;
      cout << "n, m, imax : "<<n<<" "<<m<<" "<<kmax<<endl<<endl;
      A.resize(n,m);
      vector<Triplet<double>> triplets;
      double a(0);
      int i(0);
      int j(0);
      for (int k=0; k<kmax; k++)
	{

	  flux_A >> i >> j >> a ;

    if (i!=j)
    {
	     triplets.push_back({i-1,j-1,a});
	     triplets.push_back({j-1,i-1,a});
     }

     else
     {
       triplets.push_back({i-1,j-1,a});
     }
	}
      A.setFromTriplets(triplets.begin(), triplets.end());
    }
  else
    {
      cout << "Ce fichier n'existe pas"<<endl;
    }
  flux_A.close();
  return A;
}
