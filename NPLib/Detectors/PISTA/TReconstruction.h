#ifndef TReconstruction_h
#define TReconstruction_h

using namespace std;

class TReconstruction{
  public:
    TReconstruction();
    ~TReconstruction();

  public:
    void ReadMatrix(Char_t* FName);

  private:
    float GetRandom();
    void RandomInit();
      
  private: 

    ClassDef(TReconstruction,1)
};

#endif

