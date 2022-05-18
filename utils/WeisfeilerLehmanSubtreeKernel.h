#ifndef wl_st_kern
#define wl_st_kern

#include "md5.h"
#include <vector>
#include <map>
#include <list>
#include <string>
#include <set>
#include <algorithm>
#include <chrono>
#include <stdio.h>
#include <string.h>

typedef std::string hash_;
typedef std::pair<hash_,hash_> MNode;
typedef std::pair< MNode*,int > WeightEdge;

typedef std::vector< std::vector< WeightEdge > > WeightAdj;

//
class WLSubTreeRps{

  private:

    // Attributes
    WeightAdj MAMGraph;
    std::vector<MNode*> Nodes;
    MD5 hash_fun;
    std::map<hash_,int>::key_compare Comparator;
    std::map<hash_,int> Counter;
    long long int CountNorm;

    void readMAM(char *EdgeListIn, char *FeatsIn);
    void do_rec_even();
    void do_rec_odd();
    void updateCounter(hash_& feat);
    void load(char *CounterIn);

  public:

    WLSubTreeRps() {};
    WLSubTreeRps(
          char *EdgeListIn,
          const int depth,
          char *FeatsIn = NULL
                );

    WLSubTreeRps( char *CountIn );

    ~WLSubTreeRps();

    std::vector<double> vect_similarity(
        const std::vector<WLSubTreeRps*>& LOtherWLST,
        bool normalised = true
      );

    double similarity(
        WLSubTreeRps* LOtherWLST,
        bool normalised = true
      );

    void display();
    void display_simple();

    void save(char *CounterOut);

};


#endif //wl_st_kern
