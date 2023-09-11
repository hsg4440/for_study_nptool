#include<iostream>
#include<math.h>
#include<vector>
#include<fstream>
#include<cstdlib>
#include<ctime>
#include<string.h>
#include"recuit.h"
#include"mini.h"
using namespace std;

/////////// TEST MH 1D avec fct simple

/* int main()
{
  double ini,final,pa,lambd;
  final=1;
  ini=100;
  pa=0.1;
  lambd=0.99;
  dim1 objet(lambd,ini,final,pa);
  //cout << objet.temp1() << endl;
  objet.MetroHast();
  cout <<  endl << "Le minimum global de la fonction est situé en: "  <<  objet.read() << endl <<  endl;
  
  


  return 0;
  }*/





 ///////// TEST classe DATA fct lin



 /* int main()
{
  double a,b;
  int n;
  string fich;
  fich = "fctlin.dat";
  a= 5;
  b= 2;
  n=20;
  lin objet(a, b, n);
  objet.data::initialise(fich);
  objet.initialise();
  return 0;
  }*/

 


  ////////// TEST MH 2D pour retrouver ax + b
 

  /* int main()
{
  double a,b,ampl,etype;
  int n;
  string fich;
  fich = "fctlin.dat";
  ampl=3;
  etype=3;
  a= 5;
  b= 2;
  n=20;
  ling objet1(a, b, n,ampl,etype);
  objet1.data::initialise(fich);
  objet1.initialise();









  
  double ini,final,pa,lambd;
  final=0.0001;
  ini=100;
  pa=0.1;
  lambd=0.999;
  dim2 objet2(lambd,ini,final,pa);
  //cout << objet2.temp1() << endl;
  objet2.recuit::MetroHast(fich);

  objet2.MetroHast();
  cout << endl << "La meilleure approximation en fonction linéaire des données fournies est de la forme: " << objet2.reada() << "*x + "  << objet2.readb() <<endl << endl;
  //objet2.plot();
  
  

  return 0;} */
  

   ////////////// TEST Fonction KHI2


   /* int main()
{
  double a,b;
  int n;
  string fich;
  fich = "fctlin.dat";
  a= 5;
  b= 2;
  n=20;
  lin objet(a, b, n);
  objet.data::initialise(fich);
  objet.initialise();

  mcarres obj(5,2);
  obj.minimisation::initialise(fich);
  obj.initialise();
  cout << obj.read() << endl;
  return 0;
  }*/
 

    ///////// TEST ling


    /*  int main()
{
  double a,b, ampl, etype;
  int n;
  string fich;
  fich = "fctlin.dat";
  a= 5;
  b= 2;
  n=20;
  ampl= 5;
  etype = 1;
  ling objet(a, b, n, ampl,etype);
  objet.data::initialise(fich);
  objet.initialise();
  return 0;
  }
    */


     //////////// TEST ax+b avec bruit gaussien

     /*  int main()
{
  double a,b,ampl1,etype;
  int n;
  string fich;
  fich = "fctlin.dat";
  a= 5;
  b= 2;
  n=100;
  ampl1=1;
  etype=1;
  ling objet1(a, b, n, ampl1, etype);
  objet1.data::initialise(fich);
  objet1.initialise();
  return 0;
}

     */

  ///////////// TEST data Voyageur de commerce

     

      /* int main()
{
  double lmin,lmax;
  int n;
  string fich;
  fich = "fctlin.dat";
  lmin=1;
  lmax=2;
  n=4;
  tab objet2(n, lmin, lmax);
  objet2.data::initialise(fich);
  objet2.initialise();
  return 0;
}
      
      */
      


       //////// TEST dist

      
       /*  int main()
{
  double lmin,lmax;
  int n;
  string fich;
  fich = "fctlin.dat";
  lmin=1;
  lmax=2;
  n=4;
  tab objet(n, lmin, lmax);
  objet.data::initialise(fich);
  objet.initialise();

  int*li=(int*)malloc(n*sizeof(int));
  li[0]=0; li[1]=1; li[2]=3; li[3]=2;
  dist obj(n,li);
  obj.minimisation::initialise(fich);
  obj.set();
  obj.initialise();
  cout << obj.read() << endl;
  double**t=obj.readm();
  cout << t[2][0] << endl;
  return 0;
  }*/
       
  
	/////// TEST Voyageur de commerce

		 int main()
{
  double lmin,lmax;
  int n;
  string fich;
  fich = "fctlin.dat";
  lmin=1;
  lmax=2;
  n=40;
  tab objet(n, lmin, lmax);
  objet.data::initialise(fich);
  objet.initialise();
	








  
  double ini,final,lambd;
  final=0.01;
  ini=1;
  lambd=0.99;
  voy objet2(lambd,ini,final,n);
  objet2.recuit::MetroHast(fich);

  objet2.MetroHast();
  cout << "La distance minimale à parcourir est: " <<  objet2.read() << endl << endl;
  objet2.plot();
  
  


  return 0;} 
