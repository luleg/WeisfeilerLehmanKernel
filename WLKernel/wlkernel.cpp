#include "WeisfeilerLehmanSubtreeKernel.h"


char *File1,*File2;
char *Feat1 = NULL;
char *Feat2 = NULL;
int depth = 2;
bool save = false;
bool normalised = true;

void ReadArgs(int argc, char **argv){
  if (argc<3){
    fprintf(stderr,"At least the paths to 2 files containing edgelists must be provided.\n");
    exit(EXIT_FAILURE);
  }
  File1 = argv[1];
  File2 = argv[2];
  int num_arg = 3;
  while( num_arg < argc){
    if(strcmp(argv[num_arg],"-d") == 0){
      num_arg++;
      depth = atoi(argv[num_arg]);
    }
    else if(strcmp(argv[num_arg],"-save") == 0){
      save=true;
    }
    else if(strcmp(argv[num_arg],"-nonorm") == 0){
      normalised = false;
    }
    else if(strcmp(argv[num_arg],"-feat1") == 0){
      num_arg++;
      Feat1 = argv[num_arg];
    }
    else if(strcmp(argv[num_arg],"-feat2") == 0){
      num_arg++;
      Feat2 = argv[num_arg];
    }
    else{
      fprintf(stderr,"Unknown argument: %s\n",argv[num_arg]);
    }

    num_arg++;
  }
}

int main(int argc, char **argv){

  ReadArgs(argc,argv);

  WLSubTreeRps WeisLehm1,WeisLehm2;
  char *ptr = strstr(File1,".emb");
  if (ptr == NULL){
    fprintf(stderr, "Weisfeiler Lehman embedding of %s built from the graph.\n",File1);
    WeisLehm1 = WLSubTreeRps(File1,depth,Feat1);
    if (save){
      char nFile[strlen(File1)+4];
      strcpy(nFile,File1);
      strcat(nFile,".emb");
      WeisLehm1.save(nFile);
      fprintf(stderr, "Embedding saved at %s.\n",nFile);
    }
  }
  else {
    fprintf(stderr, "Weisfeiler Lehman embedding of %s directly read.\n",File1);
    WeisLehm1 = WLSubTreeRps(File1);
  }

  ptr = strstr(File2,".emb");
  if (ptr == NULL){
    fprintf(stderr, "Weisfeiler Lehman embedding of %s built from the graph.\n",File2);
    WeisLehm2 = WLSubTreeRps(File2,depth,Feat2);
    if (save){
      char nFile[strlen(File2)+4];
      strcpy(nFile,File2);
      strcat(nFile,".emb");
      WeisLehm2.save(nFile);
      fprintf(stderr, "Embedding saved at %s.\n",nFile);
    }
  }
  else {
    fprintf(stderr, "Weisfeiler Lehman embedding of %s directly read.\n",File2);
    WeisLehm2 = WLSubTreeRps(File2);
  }
  double sim = WeisLehm1.similarity(&WeisLehm2,normalised);
  fprintf(stderr,"similarity=\n");
  fprintf(stdout,"%.4e\n",sim);



  return 0;
}
