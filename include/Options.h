#ifndef DRUID_ARGS_H
#define DRUID_ARGS_H

struct Options {

    double MCPtCut;
    bool printVersion;
    bool printHelp;
    
    void parse(int &argc, char** &argv);

};

extern Options gOptions;


#endif // DRUID_ARGS_H

