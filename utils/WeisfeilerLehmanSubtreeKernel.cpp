#include "WeisfeilerLehmanSubtreeKernel.h"

// std::string hash_fun(const std::string& a_str){
//   return a_str;
// }
//
std::chrono::steady_clock::time_point start_t,end_t;

void tic(){
   start_t = std::chrono::steady_clock::now();
}


void toc(){
   end_t = std::chrono::steady_clock::now();
   std::chrono::duration<double> span = std::chrono::duration_cast<std::chrono::duration<double>>(end_t-start_t);
   fprintf(stderr,"elapsed time = %fs.\n",span);
}

void WLSubTreeRps::readMAM(char *EdgeListIn, char *FeatsIn){
  uint nb_nodes = 0;

  if (FeatsIn){
    FILE * featsinput = fopen(FeatsIn,"r");
    if (featsinput == NULL){
      fprintf(stderr,"The file %s does not exist\n", FeatsIn);
      exit(EXIT_FAILURE);
    }
    while (!feof(featsinput)) {
      uint node;
      char tmp[80];
      int nbRead = fscanf(featsinput,"%d %s\n",&node,&tmp);
      hash_ feat = hash_fun(tmp);
      hash_ lab = std::to_string(node);
      if (nbRead == 2){
        uint mx_crt = std::max(nb_nodes,node+1);
        if (mx_crt > nb_nodes){
          nb_nodes = mx_crt;
          MAMGraph.resize(nb_nodes);
          Nodes.resize(nb_nodes);
        }
        Nodes[node] = new MNode(std::make_pair(feat,lab));
        updateCounter(feat);

      }
    }
    fclose(featsinput);
  }

  MAMGraph.resize(nb_nodes);
  std::vector<int> node_degs;

  FILE * finput = fopen(EdgeListIn,"r");
  if (finput == NULL){
    fprintf(stderr,"The file %s does not exist\n", EdgeListIn);
    exit(EXIT_FAILURE);
  }
  while (!feof(finput)) {
    uint src,tgt,weight;
    int nbRead = fscanf(finput,"%d %d %d\n",&src,&tgt,&weight);
    if (src>tgt){
      continue;
    }
    if (nbRead == 3){
      uint mx_crt = std::max(src,tgt)+1;
      if (mx_crt > nb_nodes){
        nb_nodes = mx_crt;
        MAMGraph.resize(nb_nodes);
        Nodes.resize(nb_nodes);
        node_degs.resize(nb_nodes);
      }
      if (not FeatsIn){
        if (MAMGraph[src].empty()){ // Create src Node
          hash_ tmp;
          Nodes[src] = new MNode(std::make_pair(tmp,tmp));
          node_degs[src] = weight;

        }
        else{
          node_degs[src] += weight;
        }
        if (MAMGraph[tgt].empty()){ // Create tgt Node
          hash_ tmp;
          Nodes[tgt] = new MNode(std::make_pair(tmp,tmp));
          node_degs[tgt] = weight;
        }
        else{
          node_degs[tgt] += weight;
        }
      }
      WeightEdge ToTgt(Nodes[tgt],weight);
      if (MAMGraph[src].empty()){
        MAMGraph[src].resize(1);
        MAMGraph[src][0] = ToTgt;
      }
      else{
        MAMGraph[src].push_back(ToTgt);
      }
      WeightEdge ToSrc(Nodes[src],weight);
      if (MAMGraph[tgt].empty()){
        MAMGraph[tgt].resize(1);
        MAMGraph[tgt][0] = ToSrc;
      }
      else{
        MAMGraph[tgt].push_back(ToSrc);
      }
    }
  }
  fclose(finput);

  if (not FeatsIn){
    for (int node = 0; node < nb_nodes; node ++){
      if (not MAMGraph[node].empty()){
        hash_ feat = hash_fun(std::to_string(node_degs[node]));
        hash_ lab  = std::to_string(node);
        Nodes[node]->first = feat;
        Nodes[node]->second = lab;
        updateCounter(feat);
      }
    }
  }
}

WLSubTreeRps::WLSubTreeRps( char *CountIn ){
  load(CountIn);
  CountNorm = 0;
  for (std::map<hash_,int>::iterator it = Counter.begin(); it!= Counter.end(); it ++){
    CountNorm+=(it->second)*(it->second);
  }
  if (CountNorm == 0){
    CountNorm = -1;
  }
  Comparator = Counter.key_comp();
}

WLSubTreeRps::WLSubTreeRps(
                char *EdgeListIn,
                const int depth,
                char *FeatsIn
                          ) {

  readMAM(EdgeListIn,  FeatsIn);
  //
  for (int deep = 1; deep <= depth; deep++){
    fprintf(stderr, "%dth iteration:\t",deep);
    if (deep%2 == 1){
      do_rec_even();
    }
    else {
      do_rec_odd();
    }
  }
  for (int node = 0; node < Nodes.size();node ++){
    delete Nodes[node];
  }
  MAMGraph.clear();
  Nodes.clear();

  CountNorm = 0;
  for (std::map<hash_,int>::iterator it = Counter.begin(); it!= Counter.end(); it ++){
    CountNorm+=(it->second)*(it->second);
  }
  if (CountNorm == 0){
    CountNorm = -1;
  }
  Comparator = Counter.key_comp();

}
//
WLSubTreeRps::~WLSubTreeRps(){

  // for (int node = 0; node < Nodes.size();node ++){
  //   delete Nodes[node];
  // }
}
//
// This is for debbugging
void WLSubTreeRps::display(){
  fprintf(stderr, "Display the Graph with new labels\n" );
  for (int node = 0; node < Nodes.size(); node ++){
    if (not MAMGraph[node].empty()){
      fprintf(stderr,"node %d (%s,%s) : \n",node,Nodes[node]->second.c_str(),Nodes[node]->first.c_str());
      for (int nei = 0; nei < MAMGraph[node].size(); nei++){
        fprintf(stderr,"\t %s(%s) -- %d\n",MAMGraph[node][nei].first->second.c_str(),MAMGraph[node][nei].first->first.c_str(),MAMGraph[node][nei].second);
      }
    }
  }
  fprintf(stderr, "Display the multiset\n" );
  for (std::map<hash_,int>::iterator it = Counter.begin(); it!= Counter.end(); it ++){
    fprintf(stderr,"%s %d\n",it->first.c_str(),it->second);
  }
  fprintf(stderr,"And the norm? %d\n",CountNorm);
}

void WLSubTreeRps::display_simple(){
  fprintf(stderr, "Display the first 10 elements of Counter (of size) \n", Counter.size());
  std::map<hash_,int>::iterator it = Counter.begin();
  for (int u = 0; u<10; u++){
    fprintf(stderr,"%s %d\n",it->first.c_str(),it->second);
    it++;
  }
  fprintf(stderr,"And the norm? %d\n",CountNorm);
}

void WLSubTreeRps::do_rec_even(){
  tic();
  for (int node = 0; node < Nodes.size(); node++){
    if (not MAMGraph[node].empty()){
      std::vector< WeightEdge >& Adj = MAMGraph[node];
      std::multiset<hash_> multiset_lab;
      for (int nei = 0; nei < Adj.size(); nei++ ){
        std::multiset<hash_>::iterator iter_set = multiset_lab.insert(Adj[nei].first->first);
        for ( int w = 1; w< Adj[nei].second; w ++ ){
          iter_set = multiset_lab.insert(iter_set,Adj[nei].first->first);
        }
      }
      hash_ new_lab = Nodes[node]->first;
      for (std::multiset<hash_>::iterator iter = multiset_lab.begin(); iter!= multiset_lab.end();iter ++){
        new_lab+="_"+*iter;
      }
      hash_ feat = hash_fun(new_lab);
      Nodes[node]->second = feat;
      updateCounter(feat);
    }
  }
  toc();
}

void WLSubTreeRps::do_rec_odd(){
  tic();
  for (int node = 0; node < Nodes.size(); node++){
    if (not MAMGraph[node].empty()){
      std::vector< WeightEdge >& Adj = MAMGraph[node];
      std::multiset<hash_> multiset_lab;
      for (int nei = 0; nei < Adj.size(); nei++ ){
        std::multiset<hash_>::iterator iter_set = multiset_lab.insert(Adj[nei].first->second);
        for ( int w = 1; w< Adj[nei].second; w ++ ){
          iter_set = multiset_lab.insert(iter_set,Adj[nei].first->second);
        }
      }
      hash_ new_lab = Nodes[node]->second;
      for (std::multiset<hash_>::iterator iter = multiset_lab.begin(); iter!= multiset_lab.end();iter ++){
        new_lab+=" "+*iter;
      }
      hash_ feat = hash_fun(new_lab);
      Nodes[node]->first = feat;
      updateCounter(feat);
    }
  }
  toc();
}


void WLSubTreeRps::updateCounter(hash_& feat){
  std::pair<std::map<hash_,int>::iterator,bool> it_count = Counter.insert(std::pair<hash_,int>(feat,1));
  if (not it_count.second){
    it_count.first->second += 1;
  }
}

std::vector<double> WLSubTreeRps::vect_similarity(
                              const std::vector<WLSubTreeRps*>& LOtherWLST,
                              bool normalised){

  std::vector<double> res(LOtherWLST.size(),0.0);

  int max_norm = 1;

  for (int num = 0; num < LOtherWLST.size(); num ++){
    tic();
    WLSubTreeRps* crtWLST = LOtherWLST[num];
    res[num] = similarity(crtWLST,normalised);
    toc();
  }
  return res;
}

double WLSubTreeRps::similarity(
                WLSubTreeRps* OtherWLST,
                bool normalised){

  fprintf(stderr,"Compute similarity:\t");
  tic();
  int max_norm(1.0);
  if (normalised){
    max_norm = std::max(CountNorm,OtherWLST->CountNorm);
  }
  std::map<hash_,int>& Count_base = Counter;
  std::map<hash_,int>& Count_oth  = OtherWLST->Counter;
  if (Count_base.size() > Count_oth.size() ){
    Count_base = Count_oth;
    Count_oth  = Counter;
  }
  double accum = 0.0;

  std::map<hash_,int>::iterator it_oth = Count_oth.begin();
  for (std::map<hash_,int>::iterator it_base = Count_base.begin(); it_base != Count_base.end(); it_base ++){
    while(it_oth != Count_oth.end() and Comparator(it_oth->first,it_base->first)){ // this means it_oth < it_base
      it_oth ++;
    }
    if (it_oth == Count_oth.end()){
      break;
    }
    if (it_oth->first == it_base->first){
      accum += 1.0*(it_oth->second)*(it_base->second)/max_norm;
      it_oth ++;
    }
  }
  toc();
  return accum;
}

void WLSubTreeRps::save(char *CounterOut){
  FILE * foutput = fopen(CounterOut,"w");
  if (foutput == NULL){
    fprintf(stderr,"The file %s cannot be opened\n", CounterOut);
    exit(EXIT_FAILURE);
  }
  for (std::map<hash_,int>::iterator it = Counter.begin(); it!= Counter.end();it ++){
    fprintf(foutput,"%d %s\n",it->second,it->first.c_str());
  }
  fclose(foutput);
}

void WLSubTreeRps::load(char *CounterIn){
  FILE * finput = fopen(CounterIn,"r");
  if (finput == NULL){
    fprintf(stderr,"The file %s does not exist.\n", CounterIn);
    exit(EXIT_FAILURE);
  }
  std::map<hash_,int>::iterator it_emplace = Counter.begin();
  while (!feof(finput)) {
    uint occ;
    char label[80];
    int nbRead = fscanf(finput,"%d %s\n",&occ,&label);
    it_emplace = Counter.emplace_hint(it_emplace,std::string(label),occ);
  }

  fclose(finput);
}
