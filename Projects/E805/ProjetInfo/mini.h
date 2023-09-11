#ifndef minimisation_h
#define minimisation_h 1

#include<iostream>


// Classe minimisation qui définit plusieurs sous-classes contenant des fonctions d'optimisation

using namespace std;

class minimisation
{
 public:
  minimisation();
  ~minimisation();


  void set();   // Stockage d'une matrice contenant les éléments de fich dans mat
  void initialise();   // Calcul de la distance avec liste et mat
  double read(); 

 private:
  int n;
  double distance;
  double** mat;
  int*liste;
};
  


#endif
